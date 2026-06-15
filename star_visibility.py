# star_visibility.py

from __future__ import annotations

import csv
import io
import os
import sys
import traceback
import warnings
from dataclasses import dataclass
from datetime import date, datetime, timedelta
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from astropy.utils import iers
from PySide6.QtCore import QDate, Qt
from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication, QCheckBox, QComboBox, QDateEdit, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox, QPushButton, QPlainTextEdit, QSpinBox, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget


iers.conf.auto_download = False
iers.conf.auto_max_age = None
warnings.filterwarnings("ignore", message=".*NoResultsWarning.*")
warnings.filterwarnings("ignore", message=r".*flux\(V\).*", category=DeprecationWarning)


@dataclass(frozen=True)
class ObserverSite:
    name: str
    latitude_deg: float
    longitude_deg: float
    elevation_m: float
    timezone: str


DEFAULT_SITE = ObserverSite(name="La Mesa, CA approximate", latitude_deg=32.7678, longitude_deg=-117.0231, elevation_m=160.0, timezone="America/Los_Angeles")
CANONICAL_CATALOG_FIELDS = ["IDs", "Names", "Type", "RA(decimal hours)", "Dec(degrees)", "VMag", "RMax(arcmin)", "RMin(arcmin)", "PosAngle"]
CATALOG_HEADER = CANONICAL_CATALOG_FIELDS
CATALOG_DELIMITER = "|"
TARGET_CSV_FIELDS = ["target", "ra", "dec", "label"]
TEXT_ENCODINGS = ("utf-8-sig", "utf-8", "cp1252", "latin-1")
SIMBAD_METADATA_CACHE: dict[str, dict[str, str]] = {}
APP_VERSION = "1.0.0"

try:
    import certifi
    os.environ.setdefault("SSL_CERT_FILE", certifi.where())
    os.environ.setdefault("REQUESTS_CA_BUNDLE", certifi.where())
except Exception:
    pass


def runtime_folder() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent


def log_metadata_error(context: str, exc: Exception) -> None:
    try:
        log_path = runtime_folder() / "star_visibility_metadata_errors.log"
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(f"[{datetime.now().isoformat(timespec='seconds')}] {context}: {exc}\n")
            handle.write(traceback.format_exc())
            handle.write("\n")
    except Exception:
        pass


def read_text_with_fallback(path: str | Path) -> str:
    file_path = Path(path)
    data = file_path.read_bytes()
    last_error = None
    for encoding in TEXT_ENCODINGS:
        try:
            return data.decode(encoding).replace("\xa0", " ")
        except UnicodeDecodeError as exc:
            last_error = exc
    raise ValueError(f"Could not decode {file_path}: {last_error}")


def sniff_text_delimiter(text: str) -> str:
    first_line = next((line for line in text.splitlines() if line.strip()), "")
    if "|" in first_line:
        return "|"
    if "\t" in first_line:
        return "\t"
    return ","


def normalize_target_list_text(text: str) -> str:
    return text.replace("\t", ",").replace("|", ",")


def clean_simbad_value(value: object) -> str:
    if value is None or np.ma.is_masked(value):
        return ""
    text = str(value).strip()
    return "" if text in {"", "--", "nan", "None"} else text


def simbad_text(result, *candidate_names: str) -> str:
    normalized = {name.casefold().replace("_", "").replace(" ", ""): name for name in result.colnames}
    for candidate in candidate_names:
        column = normalized.get(candidate.casefold().replace("_", "").replace(" ", ""))
        if column:
            value = clean_simbad_value(result[column][0])
            if value:
                return value
    return ""


def format_magnitude(value: str) -> str:
    cleaned = (value or "").strip()
    if not cleaned:
        return ""
    try:
        return f"{float(cleaned):.2f}"
    except ValueError:
        return cleaned


def _normalize_prefixed_id(cleaned: str, prefix: str) -> str:
    upper = cleaned.upper()
    prefix_upper = prefix.upper()
    if upper == prefix_upper:
        return ""
    if upper.startswith(prefix_upper + " "):
        return cleaned
    if upper.startswith(prefix_upper) and len(cleaned) > len(prefix):
        next_char = cleaned[len(prefix)]
        if next_char.isdigit() or next_char in "+-J":
            return f"{prefix} {cleaned[len(prefix):].strip()}"
    return ""


def simbad_hip_id(result) -> str:
    ids_text = simbad_text(result, "IDS", "ID")
    main_id = simbad_text(result, "MAIN_ID", "main_id")
    identifiers = []
    for raw in [ids_text, main_id]:
        for part in raw.replace("\n", "|").split("|"):
            cleaned = " ".join(part.strip().split())
            if cleaned and cleaned not in identifiers:
                identifiers.append(cleaned)
    for prefix in ("HIP", "HD", "Gaia DR3", "Gaia DR2", "TYC", "BD", "CD", "CPD", "2MASS"):
        for identifier in identifiers:
            normalized = _normalize_prefixed_id(identifier, prefix)
            if normalized:
                return normalized
    return identifiers[0] if identifiers else ""


def simbad_vmag(result) -> str:
    value = simbad_text(result, "FLUX_V", "flux_V", "V", "VMAG", "Vmag")
    if value:
        return format_magnitude(value)
    for column in result.colnames:
        normalized = column.casefold().replace("_", "")
        if normalized in {"fluxv", "vmag"}:
            value = clean_simbad_value(result[column][0])
            if value:
                return format_magnitude(value)
    return ""


def lookup_target_by_name(target_name: str) -> dict:
    from astroquery.simbad import Simbad
    try:
        from astroquery.exceptions import NoResultsWarning
    except Exception:
        NoResultsWarning = Warning

    cleaned_name = target_name.strip()
    if not cleaned_name:
        raise ValueError("Target name cannot be blank.")

    simbad = Simbad()
    simbad.TIMEOUT = 4
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        warnings.simplefilter("ignore", NoResultsWarning)
        warnings.filterwarnings("ignore", message=".*NoResultsWarning.*")
        for field in ("ids", "V"):
            try:
                simbad.add_votable_fields(field)
            except Exception:
                pass
        result = simbad.query_object(cleaned_name)
    if result is None or len(result) == 0:
        raise ValueError(f"Could not find target name in SIMBAD/Sesame: {cleaned_name}")

    column_names = {name.lower(): name for name in result.colnames}
    if "ra" not in column_names or "dec" not in column_names:
        raise ValueError(f"SIMBAD result did not include RA/DEC columns. Available columns: {result.colnames}")

    ra_column = result[column_names["ra"]]
    dec_column = result[column_names["dec"]]
    ra = str(ra_column[0]).strip()
    dec = str(dec_column[0]).strip()
    ra_unit = str(getattr(ra_column, "unit", "") or "").lower()
    dec_unit = str(getattr(dec_column, "unit", "") or "").lower()

    if numeric_text(ra) and ("deg" in ra_unit or abs(float(ra)) > 24.0):
        coord = SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg)
    elif numeric_text(ra) and "deg" in dec_unit and ":" not in ra and "h" not in ra.lower():
        coord = SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg)
    else:
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
    degrees = coord_to_decimal_degrees(coord)
    return target_record(label=cleaned_name, target_name=cleaned_name, ra=ra, dec=dec, coord=coord, hip_id=simbad_hip_id(result), magnitude=simbad_vmag(result)) | degrees

def numeric_text(value: str) -> bool:
    try:
        float(value.strip())
        return True
    except ValueError:
        return False


def parse_target_coordinates(ra: str, dec: str, label: str | None = None) -> dict:
    ra_clean = ra.strip(); dec_clean = dec.strip()
    if not ra_clean or not dec_clean:
        raise ValueError("RA and DEC cannot be blank.")

    ra_lower = ra_clean.lower()
    if "h" in ra_lower or ":" in ra_clean:
        coord = SkyCoord(ra=ra_clean, dec=dec_clean, unit=(u.hourangle, u.deg))
    elif ra_lower.endswith("d"):
        coord = SkyCoord(ra=ra_clean, dec=dec_clean)
    elif numeric_text(ra_clean):
        coord = SkyCoord(ra=float(ra_clean) * u.deg, dec=float(dec_clean) * u.deg)
    else:
        coord = SkyCoord(ra=ra_clean, dec=dec_clean, unit=(u.hourangle, u.deg))

    return target_record(label=label, target_name=None, ra=ra_clean, dec=dec_clean, coord=coord)


def target_record(label: str | None, target_name: str | None, ra: str, dec: str, coord: SkyCoord, hip_id: str = "", magnitude: str = "") -> dict:
    return {"label": label or target_name or "Manual RA/DEC target", "target_name": target_name, "ra": ra, "dec": dec, "ra_degrees": float(coord.ra.deg), "dec_degrees": float(coord.dec.deg), "hip_id": hip_id, "magnitude": magnitude, "coord": coord}

def coord_to_decimal_degrees(coord: SkyCoord) -> dict:
    return {"ra_degrees": float(coord.ra.deg), "dec_degrees": float(coord.dec.deg)}


def resolve_target(target_name: str | None = None, ra: str | None = None, dec: str | None = None, label: str | None = None) -> dict:
    if target_name is not None and target_name.strip():
        target = lookup_target_by_name(target_name)
        return {**target, "label": label or target["label"]}
    if ra is None or dec is None:
        raise ValueError("Provide either a target name, or both RA and DEC.")
    return parse_target_coordinates(ra=ra, dec=dec, label=label)


def build_local_time_grid(observing_date: str | None, site: ObserverSite, time_step_minutes: int) -> dict:
    if time_step_minutes <= 0:
        raise ValueError("Time step must be greater than zero.")

    tz = ZoneInfo(site.timezone)
    local_date = datetime.now(tz).date() if observing_date is None else date.fromisoformat(observing_date)
    local_start = datetime(local_date.year, local_date.month, local_date.day, 12, 0, 0, tzinfo=tz)
    local_end = local_start + timedelta(days=1)
    total_minutes = int((local_end - local_start).total_seconds() // 60)
    offsets = np.arange(0, total_minutes + time_step_minutes, time_step_minutes)
    local_times = [local_start + timedelta(minutes=int(offset)) for offset in offsets]
    utc_times = [local_time.astimezone(ZoneInfo("UTC")) for local_time in local_times]
    return {"local_date": local_date, "local_times": local_times, "astropy_times": Time(utc_times)}


def contiguous_windows(local_times: list[datetime], mask: np.ndarray) -> list[dict]:
    windows = []; start = None
    for idx, is_good in enumerate(mask):
        if bool(is_good) and start is None:
            start = local_times[idx]
        elif not bool(is_good) and start is not None:
            windows.append({"start": start, "end": local_times[idx - 1]}); start = None
    if start is not None:
        windows.append({"start": start, "end": local_times[-1]})
    return windows


def format_local_time(value: datetime) -> str:
    return value.strftime("%Y-%m-%d %I:%M %p %Z")


def analyze_star_visibility(
    target_name: str | None = None,
    ra: str | None = None,
    dec: str | None = None,
    label: str | None = None,
    observing_date: str | None = None,
    site: ObserverSite = DEFAULT_SITE,
    limiting_altitude_deg: float = 25.0,
    sun_altitude_limit_deg: float = -12.0,
    time_step_minutes: int = 5,
) -> dict:
    if limiting_altitude_deg < 0.0 or limiting_altitude_deg > 90.0:
        raise ValueError("Limiting altitude must be between 0 and 90 degrees.")
    if sun_altitude_limit_deg > 0.0:
        raise ValueError("Sun altitude night limit should normally be zero or negative.")

    target = resolve_target(target_name=target_name, ra=ra, dec=dec, label=label)
    time_grid = build_local_time_grid(observing_date=observing_date, site=site, time_step_minutes=time_step_minutes)
    location = EarthLocation(lat=site.latitude_deg * u.deg, lon=site.longitude_deg * u.deg, height=site.elevation_m * u.m)
    frame = AltAz(obstime=time_grid["astropy_times"], location=location)
    target_altaz = target["coord"].transform_to(frame)
    sun_altaz = get_sun(time_grid["astropy_times"]).transform_to(frame)
    target_altitude_deg = target_altaz.alt.deg
    target_azimuth_deg = target_altaz.az.deg
    observable_mask = (sun_altaz.alt.deg <= sun_altitude_limit_deg) & (target_altitude_deg >= limiting_altitude_deg)
    windows = contiguous_windows(local_times=time_grid["local_times"], mask=observable_mask)
    max_alt_idx = int(np.argmax(target_altitude_deg))
    max_observable_alt_idx = int(np.argmax(np.where(observable_mask, target_altitude_deg, -999.0))) if np.any(observable_mask) else None

    return {
        "site": {"name": site.name, "latitude_deg": site.latitude_deg, "longitude_deg": site.longitude_deg, "elevation_m": site.elevation_m, "timezone": site.timezone},
        "target": {"label": target["label"], "name": target["target_name"], "ra": target["ra"], "dec": target["dec"], "ra_degrees": round(target["ra_degrees"], 8), "dec_degrees": round(target["dec_degrees"], 8), "hip_id": target.get("hip_id", ""), "magnitude": target.get("magnitude", "")},
        "settings": {"observing_date": time_grid["local_date"].isoformat(), "limiting_altitude_deg": limiting_altitude_deg, "sun_altitude_limit_deg": sun_altitude_limit_deg, "time_step_minutes": time_step_minutes},
        "summary": {"is_observable": bool(np.any(observable_mask)), "max_altitude_deg": round(float(target_altitude_deg[max_alt_idx]), 2), "max_altitude_azimuth_deg": round(float(target_azimuth_deg[max_alt_idx]), 2), "max_altitude_time_local": format_local_time(time_grid["local_times"][max_alt_idx]), "max_observable_altitude_deg": None if max_observable_alt_idx is None else round(float(target_altitude_deg[max_observable_alt_idx]), 2), "max_observable_altitude_time_local": None if max_observable_alt_idx is None else format_local_time(time_grid["local_times"][max_observable_alt_idx])},
        "observable_windows": [{"start_local": format_local_time(window["start"]), "end_local": format_local_time(window["end"]), "duration_minutes": int((window["end"] - window["start"]).total_seconds() // 60)} for window in windows],
    }


def _target_entry(source: str, name: str = "", ra: str = "", dec: str = "", label: str = "", raw: str = "", hip_id: str = "", magnitude: str = "") -> dict:
    return {"source": source, "name": name.strip(), "ra": ra.strip(), "dec": dec.strip(), "label": (label or name).strip(), "raw": raw.strip(), "hip_id": hip_id.strip(), "magnitude": format_magnitude(magnitude)}


def parse_target_list_text(text: str) -> list[dict]:
    entries = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        cleaned = line.strip()
        if not cleaned or cleaned.startswith("#"):
            continue
        parts = [part.strip() for part in cleaned.split(",")]
        if len(parts) == 1 and parts[0]:
            entries.append(_target_entry("name", name=parts[0], label=parts[0], raw=cleaned))
        elif len(parts) == 2 and all(parts):
            entries.append(_target_entry("radec", ra=parts[0], dec=parts[1], label=cleaned, raw=cleaned))
        elif len(parts) >= 3 and parts[1] and parts[2]:
            entries.append(_target_entry("radec", name=parts[0], ra=parts[1], dec=parts[2], label=parts[0], raw=cleaned))
        else:
            raise ValueError(f"Target list line {line_number} is not usable: {cleaned}")
    if not entries:
        raise ValueError("Target list does not contain any usable targets.")
    return entries


def load_target_csv(path: str | Path) -> list[dict]:
    csv_path = Path(path)
    name_aliases = ["target", "name", "object", "ids", "label"]
    ra_aliases = ["ra", "right ascension", "ra(decimal hours)", "ra(decimal degrees)"]
    dec_aliases = ["dec", "declination", "dec(degrees)"]
    hip_aliases = ["hip", "hip id", "hip_id", "names"]
    mag_aliases = ["vmag", "magnitude", "mag", "v mag", "v_mag"]
    text = read_text_with_fallback(csv_path)
    if not text.strip():
        raise ValueError(f"Target file is empty: {csv_path}")

    delimiter = sniff_text_delimiter(text)
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    normalized = {field.strip().casefold(): field for field in (reader.fieldnames or []) if field and field.strip()}
    name_field = next((normalized[name] for name in name_aliases if name in normalized), None)
    ra_field = next((normalized[name] for name in ra_aliases if name in normalized), None)
    dec_field = next((normalized[name] for name in dec_aliases if name in normalized), None)
    hip_field = next((normalized[name] for name in hip_aliases if name in normalized), None)
    mag_field = next((normalized[name] for name in mag_aliases if name in normalized), None)

    if not name_field and not (ra_field and dec_field):
        return parse_target_list_text(normalize_target_list_text(text))

    entries = []
    for row_number, row in enumerate(reader, start=2):
        name = (row.get(name_field, "") if name_field else "").strip()
        ra = (row.get(ra_field, "") if ra_field else "").strip()
        if ra and ra_field and ra_field.strip().casefold() == "ra(decimal hours)":
            ra = f"{float(ra) * 15.0:.8f}"
        dec = (row.get(dec_field, "") if dec_field else "").strip()
        hip_id = (row.get(hip_field, "") if hip_field else "").strip()
        magnitude = (row.get(mag_field, "") if mag_field else "").strip()
        if not name and not ra and not dec:
            continue
        if ra and dec:
            entries.append(_target_entry("radec", name=name, ra=ra, dec=dec, label=name or f"{ra}, {dec}", raw=str(row), hip_id=hip_id, magnitude=magnitude))
        elif name:
            entries.append(_target_entry("name", name=name, label=name, raw=str(row), hip_id=hip_id, magnitude=magnitude))
        else:
            raise ValueError(f"Target file row {row_number} has incomplete RA/Dec and no target name.")
    if not entries:
        raise ValueError(f"Target file contains no usable targets: {csv_path}")
    return entries



def lookup_simbad_metadata(target_name: str) -> dict[str, str]:
    key = catalog_key(target_name)
    if not key:
        return {"hip_id": "", "magnitude": ""}
    if key in SIMBAD_METADATA_CACHE:
        return SIMBAD_METADATA_CACHE[key]
    try:
        target = lookup_target_by_name(target_name)
        metadata = {"hip_id": target.get("hip_id", ""), "magnitude": target.get("magnitude", "")}
    except Exception as exc:
        log_metadata_error(f"name metadata lookup failed for {target_name}", exc)
        metadata = {"hip_id": "", "magnitude": ""}
    SIMBAD_METADATA_CACHE[key] = metadata
    return metadata


def lookup_simbad_metadata_by_coord(ra: str, dec: str, radius_arcsec: float = 30.0) -> dict[str, str]:
    key = catalog_key(f"coord:{ra},{dec}")
    if not key:
        return {"hip_id": "", "magnitude": ""}
    if key in SIMBAD_METADATA_CACHE:
        return SIMBAD_METADATA_CACHE[key]
    try:
        from astroquery.simbad import Simbad
        try:
            from astroquery.exceptions import NoResultsWarning
        except Exception:
            NoResultsWarning = Warning

        coord = parse_target_coordinates(ra=ra, dec=dec)["coord"]
        simbad = Simbad()
        simbad.TIMEOUT = 4
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            warnings.simplefilter("ignore", NoResultsWarning)
            warnings.filterwarnings("ignore", message=".*NoResultsWarning.*")
            for field in ("ids", "V"):
                try:
                    simbad.add_votable_fields(field)
                except Exception:
                    pass
            result = simbad.query_region(coord, radius=radius_arcsec * u.arcsec)
        if result is None or len(result) == 0:
            metadata = {"hip_id": "", "magnitude": ""}
        else:
            metadata = {"hip_id": simbad_hip_id(result), "magnitude": simbad_vmag(result)}
    except Exception as exc:
        log_metadata_error(f"coordinate metadata lookup failed for {ra}, {dec}", exc)
        metadata = {"hip_id": "", "magnitude": ""}
    SIMBAD_METADATA_CACHE[key] = metadata
    return metadata

def enrich_single_entry_with_simbad_metadata(entry: dict) -> tuple[dict, int]:
    updated = 0
    if entry.get("hip_id") and entry.get("magnitude"):
        return entry, updated

    metadata = {"hip_id": "", "magnitude": ""}
    name = entry.get("name") or entry.get("label")
    if name:
        metadata = lookup_simbad_metadata(name)

    if (not metadata.get("hip_id") or not metadata.get("magnitude")) and entry.get("ra") and entry.get("dec"):
        coord_metadata = lookup_simbad_metadata_by_coord(entry["ra"], entry["dec"])
        metadata = {
            "hip_id": metadata.get("hip_id") or coord_metadata.get("hip_id", ""),
            "magnitude": metadata.get("magnitude") or coord_metadata.get("magnitude", ""),
        }

    if not entry.get("hip_id") and metadata.get("hip_id"):
        entry["hip_id"] = metadata["hip_id"]
        updated += 1
    if not entry.get("magnitude") and metadata.get("magnitude"):
        entry["magnitude"] = format_magnitude(metadata["magnitude"])
        updated += 1
    return entry, updated


def enrich_entries_with_simbad_metadata(entries: list[dict], max_workers: int = 8) -> int:
    from concurrent.futures import ThreadPoolExecutor, as_completed

    pending = [entry for entry in entries if not (entry.get("hip_id") and entry.get("magnitude"))]
    if not pending:
        return 0

    updated = 0
    workers = max(1, min(max_workers, len(pending)))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(enrich_single_entry_with_simbad_metadata, entry) for entry in pending]
        for future in as_completed(futures):
            try:
                _, count = future.result()
                updated += count
            except Exception as exc:
                log_metadata_error("parallel metadata worker failed", exc)
    return updated

def parse_batch_target_line(line: str) -> dict:
    return parse_target_list_text(line)[0]


def parse_target_csv(path: str | Path) -> list[str]:
    lines = []
    for entry in load_target_csv(path):
        if entry["source"] == "name":
            lines.append(entry["name"])
        elif entry["name"]:
            lines.append(f"{entry['name']}, {entry['ra']}, {entry['dec']}")
        else:
            lines.append(f"{entry['ra']}, {entry['dec']}")
    return lines


def analyze_target_entry(entry: dict, observing_date: str | None, site: ObserverSite, limiting_altitude_deg: float = 25.0, sun_altitude_limit_deg: float = -12.0, time_step_minutes: int = 5) -> dict:
    try:
        if entry["source"] == "name":
            result = analyze_star_visibility(target_name=entry["name"], label=entry.get("label") or entry["name"], observing_date=observing_date, site=site, limiting_altitude_deg=limiting_altitude_deg, sun_altitude_limit_deg=sun_altitude_limit_deg, time_step_minutes=time_step_minutes)
        elif entry["source"] == "radec":
            result = analyze_star_visibility(ra=entry["ra"], dec=entry["dec"], label=entry.get("label") or entry.get("name") or "Manual RA/DEC target", observing_date=observing_date, site=site, limiting_altitude_deg=limiting_altitude_deg, sun_altitude_limit_deg=sun_altitude_limit_deg, time_step_minutes=time_step_minutes)
            result["target"]["hip_id"] = entry.get("hip_id", "")
            result["target"]["magnitude"] = entry.get("magnitude", "")
        else:
            raise ValueError(f"Unknown target entry source: {entry.get('source')}")
        target = result["target"]; summary = result["summary"]; windows = result["observable_windows"]
        display_name = entry.get("name") or target.get("name") or target.get("label") or entry.get("label") or entry.get("raw") or "Target"
        return {"entry": entry, "ok": True, "visible": summary["is_observable"], "result": result, "target_name": display_name, "label": display_name, "hip_id": target.get("hip_id", ""), "magnitude": target.get("magnitude", ""), "ra_degrees": target["ra_degrees"], "dec_degrees": target["dec_degrees"], "max_observable_altitude_deg": summary["max_observable_altitude_deg"], "best_time_local": summary["max_observable_altitude_time_local"], "window_count": len(windows), "first_window": "None" if not windows else f"{windows[0]['start_local']} to {windows[0]['end_local']}", "error": ""}
    except Exception as exc:
        display_name = entry.get("name") or entry.get("label") or entry.get("raw") or "Target"
        return {"entry": entry, "ok": False, "visible": False, "result": None, "target_name": display_name, "label": display_name, "hip_id": "", "magnitude": "", "ra_degrees": "", "dec_degrees": "", "max_observable_altitude_deg": None, "best_time_local": None, "window_count": 0, "first_window": "", "error": str(exc)}
def analyze_target_entries(entries: list[dict], observing_date: str | None, site: ObserverSite, limiting_altitude_deg: float = 25.0, sun_altitude_limit_deg: float = -12.0, time_step_minutes: int = 5) -> list[dict]:
    return [analyze_target_entry(entry, observing_date=observing_date, site=site, limiting_altitude_deg=limiting_altitude_deg, sun_altitude_limit_deg=sun_altitude_limit_deg, time_step_minutes=time_step_minutes) for entry in entries]


def analyze_target_list(target_lines: list[str], observing_date: str | None, site: ObserverSite, limiting_altitude_deg: float = 25.0, sun_altitude_limit_deg: float = -12.0, time_step_minutes: int = 5) -> list[dict]:
    return analyze_target_entries(parse_target_list_text("\n".join(target_lines)), observing_date=observing_date, site=site, limiting_altitude_deg=limiting_altitude_deg, sun_altitude_limit_deg=sun_altitude_limit_deg, time_step_minutes=time_step_minutes)


def row_text(row: dict, *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None:
            return str(value).strip()
    return ""


def normalize_catalog_row(row: dict) -> dict:
    ids = row_text(row, "IDs", "HIP", "HIP ID")
    names = row_text(row, "Names", "Name") or ids
    object_type = row_text(row, "Type", "Object Type") or "Star"
    ra_hours = row_text(row, "RA(decimal hours)")
    if not ra_hours:
        ra_degrees = row_text(row, "RA", "RA(degrees)", "RA(decimal degrees)")
        if ra_degrees:
            ra_hours = f"{float(ra_degrees) / 15.0:.8f}"
    return {
        "IDs": ids or names,
        "Names": names,
        "Type": object_type,
        "RA(decimal hours)": ra_hours,
        "Dec(degrees)": row_text(row, "Dec(degrees)", "Dec", "DEC"),
        "VMag": format_magnitude(row_text(row, "VMag", "Magnitude")),
        "RMax(arcmin)": row_text(row, "RMax(arcmin)", "Major Axis"),
        "RMin(arcmin)": row_text(row, "RMin(arcmin)", "Minor Axis"),
        "PosAngle": row_text(row, "PosAngle", "Orientation"),
    }
def read_catalog(path: str | Path) -> list[dict]:
    catalog_path = Path(path)
    text = read_text_with_fallback(catalog_path)
    if not text.strip():
        return []
    delimiter = sniff_text_delimiter(text)
    rows = [normalize_catalog_row(row) for row in csv.DictReader(io.StringIO(text), delimiter=delimiter)]
    return [row for row in rows if any(value.strip() for value in row.values())]


def write_catalog(path: str | Path, rows: list[dict]) -> None:
    catalog_path = Path(path)
    with catalog_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CANONICAL_CATALOG_FIELDS, delimiter=CATALOG_DELIMITER, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in CANONICAL_CATALOG_FIELDS})


def append_catalog_rows(path: str | Path, rows: list[dict]) -> None:
    if not rows:
        return
    catalog_path = Path(path)
    header = CATALOG_DELIMITER.join(CATALOG_HEADER)
    existing_text = read_text_with_fallback(catalog_path) if catalog_path.exists() else ""
    if existing_text.strip():
        first_line = next((line.strip() for line in existing_text.splitlines() if line.strip()), "")
        if first_line != header:
            raise ValueError(f"Catalog first row must be exactly:\n{header}")
        needs_leading_newline = not existing_text.endswith(("\n", "\r"))
    else:
        needs_leading_newline = False

    with catalog_path.open("a", encoding="utf-8", newline="") as handle:
        if not existing_text.strip():
            handle.write(header + "\n")
        elif needs_leading_newline:
            handle.write("\n")
        writer = csv.DictWriter(handle, fieldnames=CANONICAL_CATALOG_FIELDS, delimiter=CATALOG_DELIMITER, extrasaction="ignore", lineterminator="\n")
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in CANONICAL_CATALOG_FIELDS})


def catalog_key(value: str) -> str:
    return " ".join(value.strip().casefold().split())


class StarVisibilityWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"Star Visibility Planner - CSV Target List UI v{APP_VERSION}")
        self.resize(1220, 820)
        self.target_list_text = ""
        self.target_entries: list[dict] = []
        self.target_list_source = ""
        self.target_csv_path: Path | None = None
        self.catalog_path: Path | None = None
        self.catalog_rows: list[dict] = []
        self.catalog_header: list[str] = CATALOG_HEADER.copy()
        self.last_batch_results: list[dict] = []
        self.last_metadata_summary: dict[str, int] = {"targets": 0, "hip_id": 0, "magnitude": 0, "updated": 0}

        self.target_mode = QComboBox(); self.target_mode.addItems(["Lookup by target name", "Manual RA/DEC", "Target list"])
        self.target_name = QLineEdit(); self.target_name.setPlaceholderText("Example: Vega, Betelgeuse, Sirius, M31")
        self.ra = QLineEdit(); self.ra.setPlaceholderText('Decimal degrees preferred. Also accepts "18h36m56.3s", "18:36:56.3", or "279.2347d"')
        self.dec = QLineEdit(); self.dec.setPlaceholderText('Decimal degrees, for example "+38.7837"')
        self.target_list_status = QLabel("No target CSV loaded. Use Load Target CSV... or Check Visibility to select one.")
        self.target_list_status.setWordWrap(False)
        self.target_list_status.setMinimumHeight(24)

        self.use_today = QCheckBox("Use today"); self.use_today.setChecked(True)
        self.date = QDateEdit(); self.date.setCalendarPopup(True); self.date.setDisplayFormat("yyyy-MM-dd"); self.date.setDate(QDate.currentDate()); self.date.setEnabled(False)
        self.limiting_altitude = QDoubleSpinBox(); self.limiting_altitude.setRange(0.0, 90.0); self.limiting_altitude.setDecimals(1); self.limiting_altitude.setSingleStep(1.0); self.limiting_altitude.setValue(25.0); self.limiting_altitude.setSuffix(" deg")
        self.sun_altitude = QDoubleSpinBox(); self.sun_altitude.setRange(-30.0, 0.0); self.sun_altitude.setDecimals(1); self.sun_altitude.setSingleStep(1.0); self.sun_altitude.setValue(-12.0); self.sun_altitude.setSuffix(" deg")
        self.step_minutes = QSpinBox(); self.step_minutes.setRange(1, 60); self.step_minutes.setValue(5); self.step_minutes.setSuffix(" min")

        self.site_name = QLineEdit(DEFAULT_SITE.name)
        self.latitude = QDoubleSpinBox(); self.latitude.setRange(-90.0, 90.0); self.latitude.setDecimals(6); self.latitude.setValue(DEFAULT_SITE.latitude_deg); self.latitude.setSuffix(" deg")
        self.longitude = QDoubleSpinBox(); self.longitude.setRange(-180.0, 180.0); self.longitude.setDecimals(6); self.longitude.setValue(DEFAULT_SITE.longitude_deg); self.longitude.setSuffix(" deg")
        self.elevation = QDoubleSpinBox(); self.elevation.setRange(-500.0, 9000.0); self.elevation.setDecimals(1); self.elevation.setValue(DEFAULT_SITE.elevation_m); self.elevation.setSuffix(" m")
        self.timezone = QLineEdit(DEFAULT_SITE.timezone)

        self.run_button = QPushButton("Check Visibility")
        self.load_target_csv_button = QPushButton("Load Target CSV...")
        self.load_catalog_button = QPushButton("Load Catalog CSV...")
        self.add_visible_button = QPushButton("Add Visible to Catalog"); self.add_visible_button.setEnabled(False)
        self.save_results_button = QPushButton("Save Batch CSV"); self.save_results_button.setEnabled(False)
        self.clear_button = QPushButton("Clear Results")
        self.help_button = QPushButton("Help")

        self.batch_results = QTableWidget(); self.batch_results.setColumnCount(10); self.batch_results.setHorizontalHeaderLabels(["Target", "HIP/ID", "VMag", "Visible?", "RA deg", "Dec deg", "Max Obs Alt", "Best Time", "Windows", "Error"]); self.batch_results.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch); self.batch_results.setMinimumHeight(230)
        self.output = QPlainTextEdit(); self.output.setReadOnly(True)

        self._build_layout(); self._wire_signals(); self._sync_target_mode()

    def _build_layout(self) -> None:
        target_box = QGroupBox("Target")
        target_form = QFormLayout(target_box)
        target_form.addRow("Input mode", self.target_mode)
        target_form.addRow("Target name", self.target_name)
        target_form.addRow("RA", self.ra)
        target_form.addRow("DEC", self.dec)
        target_form.addRow("List status", self.target_list_status)

        night_box = QGroupBox("Observing Settings")
        night_form = QFormLayout(night_box)
        night_form.addRow("", self.use_today)
        night_form.addRow("Observing date", self.date)
        night_form.addRow("Limiting altitude", self.limiting_altitude)
        night_form.addRow("Sun altitude night limit", self.sun_altitude)
        night_form.addRow("Time step", self.step_minutes)

        site_box = QGroupBox("Observer Site")
        site_form = QFormLayout(site_box)
        site_form.addRow("Site name", self.site_name)
        site_form.addRow("Latitude", self.latitude)
        site_form.addRow("Longitude", self.longitude)
        site_form.addRow("Elevation", self.elevation)
        site_form.addRow("Timezone", self.timezone)

        button_row = QHBoxLayout()
        for button in [self.run_button, self.load_target_csv_button, self.load_catalog_button, self.add_visible_button, self.save_results_button, self.clear_button, self.help_button]:
            button_row.addWidget(button)
        button_row.addStretch(1)

        top_row = QHBoxLayout(); top_row.addWidget(target_box); top_row.addWidget(night_box); top_row.addWidget(site_box)
        root = QVBoxLayout(); root.addLayout(top_row); root.addLayout(button_row); root.addWidget(QLabel("Batch Results")); root.addWidget(self.batch_results); root.addWidget(QLabel("Detailed Results")); root.addWidget(self.output)
        widget = QWidget(); widget.setLayout(root); self.setCentralWidget(widget)

    def _wire_signals(self) -> None:
        self.target_mode.currentIndexChanged.connect(self._sync_target_mode)
        self.use_today.toggled.connect(self._sync_date_mode)
        self.run_button.clicked.connect(self._run_analysis)
        self.load_target_csv_button.clicked.connect(self._load_target_csv)
        self.load_catalog_button.clicked.connect(self._load_catalog)
        self.add_visible_button.clicked.connect(self._add_visible_to_catalog)
        self.save_results_button.clicked.connect(self._save_batch_results_csv)
        self.clear_button.clicked.connect(self._clear_results)
        self.help_button.clicked.connect(self._show_help)

    def _show_help(self) -> None:
        QMessageBox.information(self, "Star Visibility Help", """Star Visibility Planner

Target list mode
- Click Load Target CSV... and select a .csv or .txt star list.
- Supported columns include target/name/object, RA, Dec, HIP/HIP ID/Names, and VMag/Magnitude/Mag.
- If RA and Dec are present, the app uses those coordinates directly.
- If only a target name is present, the app resolves the target through SIMBAD.

Identifier and magnitude
- HIP is used when available.
- If HIP is not available, the app tries HD, Gaia DR3, Gaia DR2, TYC, BD, CD, CPD, then 2MASS.
- VMag is displayed and saved to two decimal places when numeric.

Catalog append
- Catalog output is compatible with SharpCap custom catalogs.
- Load your custom catalog with Load Catalog CSV...
- Click Add Visible to Catalog after checking visibility.
- Only visible targets are appended.
- The catalog header must exactly match:
IDs|Names|Type|RA(decimal hours)|Dec(degrees)|VMag|RMax(arcmin)|RMin(arcmin)|PosAngle

Catalog output fields
- IDs = target name from the list.
- Names = HIP/alternate ID.
- RA is saved in decimal hours.
- Dec is saved in decimal degrees.
- VMag is saved immediately after Dec.""")
    def _sync_target_mode(self) -> None:
        mode = self.target_mode.currentText(); lookup_mode = mode == "Lookup by target name"; manual_mode = mode == "Manual RA/DEC"; batch_mode = mode == "Target list"
        self.target_name.setEnabled(lookup_mode); self.ra.setEnabled(manual_mode); self.dec.setEnabled(manual_mode); self.load_target_csv_button.setEnabled(batch_mode)

    def _sync_date_mode(self) -> None:
        self.date.setEnabled(not self.use_today.isChecked())

    def _make_site(self) -> ObserverSite:
        site_name = self.site_name.text().strip(); timezone = self.timezone.text().strip()
        if not site_name:
            raise ValueError("Site name cannot be blank.")
        if not timezone:
            raise ValueError("Timezone cannot be blank. Example: America/Los_Angeles")
        return ObserverSite(name=site_name, latitude_deg=float(self.latitude.value()), longitude_deg=float(self.longitude.value()), elevation_m=float(self.elevation.value()), timezone=timezone)

    def _observing_date(self) -> str | None:
        return None if self.use_today.isChecked() else self.date.date().toString("yyyy-MM-dd")

    def _get_single_inputs(self) -> dict:
        mode = self.target_mode.currentText()
        if mode == "Target list":
            raise ValueError("Target list mode should be run through batch analysis.")
        if mode == "Lookup by target name":
            target_name = self.target_name.text().strip()
            if not target_name:
                raise ValueError("Enter a target name, such as Vega, Betelgeuse, Sirius, or M31.")
            return {"target_name": target_name, "ra": None, "dec": None, "label": target_name, "observing_date": self._observing_date()}
        ra = self.ra.text().strip(); dec = self.dec.text().strip()
        if not ra or not dec:
            raise ValueError("Enter both RA and DEC, or switch to target-name lookup.")
        return {"target_name": None, "ra": ra, "dec": dec, "label": "Manual RA/DEC target", "observing_date": self._observing_date()}

    def _update_target_list_status(self) -> None:
        count = len(self.target_entries)
        if not count:
            self.target_list_status.setText("No target CSV loaded. Use Load Target CSV... or Check Visibility to select one.")
        elif self.target_list_source == "csv":
            self.target_list_status.setText(f"CSV loaded - {count} target(s) ready.")
        else:
            self.target_list_status.setText(f"{count} target(s) ready.")

    def _load_target_csv(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select Star List", "", "Star List Files (*.csv *.txt);;CSV Files (*.csv);;Text Files (*.txt);;All Files (*)")
        if not path:
            return
        try:
            entries = load_target_csv(path)
            self.target_csv_path = Path(path)
            self.target_entries = entries
            self.target_list_source = "csv"
            self.target_list_text = "\n".join(entry["name"] if entry["source"] == "name" else f"{entry.get('label') or entry.get('name') or ''}, {entry['ra']}, {entry['dec']}".strip(", ") for entry in entries)
            self.target_mode.setCurrentText("Target list")
            self._update_target_list_status()
            QMessageBox.information(self, "Target List Loaded", f"Loaded {len(entries)} target(s) from:\n{path}")
        except Exception as exc:
            QMessageBox.critical(self, "Target List Error", str(exc))

    def _load_catalog(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Load Custom Catalog", "", "Catalog/CSV Files (*.csv *.txt);;All Files (*)")
        if not path:
            return
        self.catalog_path = Path(path)
        self.catalog_rows = read_catalog(self.catalog_path)
        self.add_visible_button.setEnabled(bool(self.last_batch_results))
        QMessageBox.information(self, "Catalog Loaded", f"Loaded {len(self.catalog_rows)} catalog row(s).\n\nCatalog will be saved in the exact template format: IDs, Names, Type, RA(decimal hours), Dec(degrees), VMag, RMax(arcmin), RMin(arcmin), PosAngle.")

    def _run_analysis(self) -> None:
        try:
            if self.target_mode.currentText() == "Target list":
                self._run_batch_analysis(); return
            inputs = self._get_single_inputs()
            result = analyze_star_visibility(target_name=inputs["target_name"], ra=inputs["ra"], dec=inputs["dec"], label=inputs["label"], observing_date=inputs["observing_date"], site=self._make_site(), limiting_altitude_deg=float(self.limiting_altitude.value()), sun_altitude_limit_deg=float(self.sun_altitude.value()), time_step_minutes=int(self.step_minutes.value()))
            self.last_batch_results = []
            self.batch_results.setRowCount(0)
            self.save_results_button.setEnabled(False)
            self.add_visible_button.setEnabled(False)
            self.output.setPlainText(self._format_result(result))
        except Exception as exc:
            QMessageBox.critical(self, "Visibility Analysis Error", f"{exc}\n\nDetails:\n{traceback.format_exc()}")

    def _run_batch_analysis(self) -> None:
        if not self.target_entries:
            self._load_target_csv()
            if not self.target_entries:
                raise ValueError("Select a target CSV before running Target list mode.")
        self.target_list_status.setText("Resolving ID/VMag metadata from SIMBAD...")
        QApplication.processEvents()
        updated_metadata = enrich_entries_with_simbad_metadata(self.target_entries)
        self.last_metadata_summary = {"targets": len(self.target_entries), "hip_id": sum(1 for entry in self.target_entries if entry.get("hip_id")), "magnitude": sum(1 for entry in self.target_entries if entry.get("magnitude")), "updated": updated_metadata}
        self.target_list_status.setText(f"CSV loaded - {len(self.target_entries)} target(s). Metadata: HIP/ID {self.last_metadata_summary['hip_id']}/{len(self.target_entries)}, VMag {self.last_metadata_summary['magnitude']}/{len(self.target_entries)}.")
        QApplication.processEvents()
        results = analyze_target_entries(entries=self.target_entries, observing_date=self._observing_date(), site=self._make_site(), limiting_altitude_deg=float(self.limiting_altitude.value()), sun_altitude_limit_deg=float(self.sun_altitude.value()), time_step_minutes=int(self.step_minutes.value()))
        self.last_batch_results = results
        self._populate_batch_results(results)
        self.output.setPlainText(self._format_batch_report(results))
        self.save_results_button.setEnabled(bool(results))
        self.add_visible_button.setEnabled(bool(results))

    def _populate_batch_results(self, results: list[dict]) -> None:
        self.batch_results.setRowCount(len(results))
        for row_idx, item in enumerate(results):
            values = [item["target_name"], item.get("hip_id", "") or "--", item.get("magnitude", "") or "--", "Yes" if item["visible"] else "No", str(item["ra_degrees"]), str(item["dec_degrees"]), "--" if item["max_observable_altitude_deg"] is None else f"{item['max_observable_altitude_deg']} deg", "--" if item["best_time_local"] is None else str(item["best_time_local"]), str(item["window_count"]) if item["ok"] else "--", item["error"]]
            for col_idx, value in enumerate(values):
                table_item = QTableWidgetItem(value); table_item.setFlags(table_item.flags() & ~Qt.ItemFlag.ItemIsEditable); self.batch_results.setItem(row_idx, col_idx, table_item)

    def _format_result(self, result: dict) -> str:
        site = result["site"]; target = result["target"]; settings = result["settings"]; summary = result["summary"]
        lines = ["Star Visibility Report", "======================", "", f"Site: {site['name']}", f"Location: lat {site['latitude_deg']:.6f} deg, lon {site['longitude_deg']:.6f} deg, elevation {site['elevation_m']:.1f} m", f"Timezone: {site['timezone']}", f"Observing date: {settings['observing_date']}", "", f"Target: {target['label']}", f"HIP/ID: {target.get('hip_id', '') or '--'}", f"VMag: {target.get('magnitude', '') or '--'}", f"RA input: {target['ra']}", f"DEC input: {target['dec']}", f"RA decimal degrees: {target['ra_degrees']:.8f}", f"Dec decimal degrees: {target['dec_degrees']:.8f}", "", f"Limiting altitude: {settings['limiting_altitude_deg']} deg", f"Night definition: Sun altitude <= {settings['sun_altitude_limit_deg']} deg", f"Time step: {settings['time_step_minutes']} min", "", f"Maximum altitude over noon-to-noon interval: {summary['max_altitude_deg']} deg", f"Azimuth at maximum altitude: {summary['max_altitude_azimuth_deg']} deg", f"Time of maximum altitude: {summary['max_altitude_time_local']}", ""]
        if not summary["is_observable"]:
            lines.append("Result: NOT observable at night above the limiting altitude."); return "\n".join(lines)
        lines.extend(["Result: Observable.", f"Maximum observable altitude: {summary['max_observable_altitude_deg']} deg at {summary['max_observable_altitude_time_local']}", "", "Observable local-time windows:"])
        for window in result["observable_windows"]:
            lines.append(f"  {window['start_local']} to {window['end_local']}  ({window['duration_minutes']} min)")
        return "\n".join(lines)

    def _format_batch_report(self, results: list[dict]) -> str:
        metadata = self.last_metadata_summary
        lines = ["Batch Star Visibility Report", "============================", "", f"Targets checked: {len(results)}", f"Visible targets: {sum(1 for item in results if item['visible'])}", f"Lookup/errors: {sum(1 for item in results if not item['ok'])}", f"Metadata filled: HIP/ID {metadata.get('hip_id', 0)}/{metadata.get('targets', len(results))}, VMag {metadata.get('magnitude', 0)}/{metadata.get('targets', len(results))}", f"Metadata fields updated this run: {metadata.get('updated', 0)}", "", "Each target below includes the full discrete visibility details written to the output.", "Output units: RA decimal degrees, Dec decimal degrees.", ""]
        for index, item in enumerate(results, start=1):
            entry = item.get("entry", {})
            lines.extend(["-" * 72, f"Target {index} of {len(results)}: {item['target_name']}", f"List source: {entry.get('source', '')}"])
            if entry.get("raw"):
                lines.append(f"List row: {entry['raw']}")
            if entry.get("name"):
                lines.append(f"Input name: {entry['name']}")
            if entry.get("ra") or entry.get("dec"):
                lines.append(f"Input RA/Dec: {entry.get('ra', '')} / {entry.get('dec', '')}")
            lines.append("")
            if not item["ok"]:
                lines.extend(["Result: ERROR", f"Error detail: {item['error']}", ""])
                continue
            lines.append(self._format_result(item["result"]))
            lines.append("")
        return "\n".join(lines).rstrip() + "\n"

    def _batch_csv_row(self, item: dict) -> dict:
        entry = item.get("entry", {})
        result = item.get("result") or {}
        summary = result.get("summary", {})
        target = result.get("target", {})
        windows = result.get("observable_windows", [])
        all_windows = "; ".join(f"{window['start_local']} to {window['end_local']} ({window['duration_minutes']} min)" for window in windows)
        return {
            "target_name": item.get("target_name", ""),
            "source": entry.get("source", ""),
            "input_name": entry.get("name", ""),
            "input_ra": entry.get("ra", ""),
            "input_dec": entry.get("dec", ""),
            "raw_input": entry.get("raw", ""),
            "ok": item.get("ok", ""),
            "visible": item.get("visible", ""),
            "resolved_label": target.get("label", ""),
            "resolved_name": target.get("name", ""),
            "resolved_ra_input": target.get("ra", ""),
            "resolved_dec_input": target.get("dec", ""),
            "hip_id": target.get("hip_id", item.get("hip_id", "")),
            "magnitude": target.get("magnitude", item.get("magnitude", "")),
            "ra_degrees": item.get("ra_degrees", ""),
            "dec_degrees": item.get("dec_degrees", ""),
            "max_altitude_deg": summary.get("max_altitude_deg", ""),
            "max_altitude_azimuth_deg": summary.get("max_altitude_azimuth_deg", ""),
            "max_altitude_time_local": summary.get("max_altitude_time_local", ""),
            "max_observable_altitude_deg": item.get("max_observable_altitude_deg", ""),
            "best_time_local": item.get("best_time_local", ""),
            "first_window": item.get("first_window", ""),
            "all_windows": all_windows,
            "window_count": item.get("window_count", ""),
            "error": item.get("error", ""),
        }

    def _save_batch_results_csv(self) -> None:
        if not self.last_batch_results:
            raise ValueError("No batch results to save.")
        path, _ = QFileDialog.getSaveFileName(self, "Save Batch Results CSV", "star_visibility_results.csv", "CSV Files (*.csv);;All Files (*)")
        if not path:
            return
        fieldnames = ["target_name", "source", "input_name", "input_ra", "input_dec", "raw_input", "ok", "visible", "resolved_label", "resolved_name", "resolved_ra_input", "resolved_dec_input", "hip_id", "magnitude", "ra_degrees", "dec_degrees", "max_altitude_deg", "max_altitude_azimuth_deg", "max_altitude_time_local", "max_observable_altitude_deg", "best_time_local", "first_window", "all_windows", "window_count", "error"]
        with Path(path).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for item in self.last_batch_results:
                writer.writerow(self._batch_csv_row(item))
        QMessageBox.information(self, "Saved", f"Saved batch results to:\n{path}")
    def _add_visible_to_catalog(self) -> None:
        if self.catalog_path is None:
            self._load_catalog()
            if self.catalog_path is None:
                return
        self.catalog_rows = read_catalog(self.catalog_path)
        existing = {catalog_key(row.get("IDs", "") or row.get("Names", "")) for row in self.catalog_rows}
        new_rows = []
        skipped = 0
        for item in self.last_batch_results:
            if not item["ok"] or not item["visible"]:
                continue
            display_name = item["target_name"]
            name_key = catalog_key(display_name)
            if name_key in existing:
                skipped += 1
                continue
            target = item.get("result", {}).get("target", {})
            hip_id = target.get("hip_id", "") or item.get("hip_id", "")
            magnitude = format_magnitude(target.get("magnitude", "") or item.get("magnitude", ""))
            row = {"IDs": display_name, "Names": hip_id, "Type": "Star", "RA(decimal hours)": f"{float(item['ra_degrees']) / 15.0:.8f}", "Dec(degrees)": f"{float(item['dec_degrees']):.6f}", "VMag": magnitude, "RMax(arcmin)": "", "RMin(arcmin)": "", "PosAngle": ""}
            new_rows.append(row)
            self.catalog_rows.append(row)
            existing.add(name_key)
        append_catalog_rows(self.catalog_path, new_rows)
        QMessageBox.information(self, "Catalog Updated", f"Added {len(new_rows)} visible target(s) to catalog.\nSkipped {skipped} duplicate(s).\nCatalog appended:\n{self.catalog_path}")

    def _clear_results(self) -> None:
        self.output.clear(); self.batch_results.setRowCount(0); self.last_batch_results = []; self.save_results_button.setEnabled(False); self.add_visible_button.setEnabled(False)


def apply_night_vision_theme(app: QApplication) -> None:
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor(5, 0, 0)); palette.setColor(QPalette.ColorRole.WindowText, QColor(255, 80, 80)); palette.setColor(QPalette.ColorRole.Base, QColor(0, 0, 0)); palette.setColor(QPalette.ColorRole.AlternateBase, QColor(20, 0, 0)); palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(20, 0, 0)); palette.setColor(QPalette.ColorRole.ToolTipText, QColor(255, 100, 100)); palette.setColor(QPalette.ColorRole.Text, QColor(255, 80, 80)); palette.setColor(QPalette.ColorRole.Button, QColor(25, 0, 0)); palette.setColor(QPalette.ColorRole.ButtonText, QColor(255, 90, 90)); palette.setColor(QPalette.ColorRole.BrightText, QColor(255, 0, 0)); palette.setColor(QPalette.ColorRole.Highlight, QColor(120, 0, 0)); palette.setColor(QPalette.ColorRole.HighlightedText, QColor(255, 180, 180)); app.setPalette(palette)
    app.setStyleSheet("""
        QWidget { background-color: #050000; color: #ff5050; selection-background-color: #780000; selection-color: #ffc0c0; font-size: 10pt; }
        QGroupBox { border: 1px solid #7a1a1a; border-radius: 4px; margin-top: 10px; padding-top: 10px; color: #ff6666; }
        QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 4px 0 4px; color: #ff6666; }
        QLineEdit, QPlainTextEdit, QComboBox, QDateEdit, QDoubleSpinBox, QSpinBox, QTableWidget { background-color: #000000; color: #ff5555; border: 1px solid #772020; border-radius: 4px; padding: 4px; }
        QHeaderView::section { background-color: #240000; color: #ff6666; border: 1px solid #772020; padding: 4px; }
        QLineEdit:disabled, QDateEdit:disabled { background-color: #100000; color: #804040; border: 1px solid #401010; }
        QPushButton { background-color: #240000; color: #ff6666; border: 1px solid #993333; border-radius: 4px; padding: 6px 12px; }
        QPushButton:hover { background-color: #3a0000; border: 1px solid #cc4444; }
        QPushButton:pressed { background-color: #5a0000; }
        QPushButton:disabled { background-color: #120000; color: #704040; border: 1px solid #401010; }
        QCheckBox { color: #ff6666; }
        QCheckBox::indicator { width: 16px; height: 16px; }
        QCheckBox::indicator:unchecked { background-color: #000000; border: 1px solid #993333; }
        QCheckBox::indicator:checked { background-color: #aa0000; border: 1px solid #ff5555; }
        QMessageBox { background-color: #050000; color: #ff5555; }
    """)


def main() -> None:
    app = QApplication.instance()
    owns_app = app is None
    if app is None:
        app = QApplication(sys.argv)
    apply_night_vision_theme(app)
    window = StarVisibilityWindow()
    window.show(); window.raise_(); window.activateWindow()
    app._star_visibility_window = window
    if owns_app:
        sys.exit(app.exec())


if __name__ == "__main__":
    main()













