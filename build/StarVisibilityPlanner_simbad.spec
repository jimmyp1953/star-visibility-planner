# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_data_files
from PyInstaller.utils.hooks import copy_metadata
from pathlib import Path
import sys

datas = []
datas += collect_data_files('astropy')
datas += collect_data_files('astropy_iers_data')
datas += collect_data_files('astroquery')
datas += copy_metadata('astroquery')

project_root = Path(SPECPATH).parent.parent
env_root = Path(sys.prefix)
conda_bin = env_root / 'Library' / 'bin'
conda_dlls = env_root / 'DLLs'
qt_dll_names = [
    'pyside6.cp313-win_amd64.dll',
    'shiboken6.cp313-win_amd64.dll',
    'Qt6Core.dll',
    'Qt6Gui.dll',
    'Qt6Widgets.dll',
    'Qt6Network.dll',
    'Qt6Svg.dll',
    'libssl-3-x64.dll',
    'libcrypto-3-x64.dll',
]
binaries = [(str(conda_bin / name), '.') for name in qt_dll_names if (conda_bin / name).exists()]
binaries += [(str(conda_dlls / name), '.') for name in ['_ssl.pyd', '_hashlib.pyd'] if (conda_dlls / name).exists()]


a = Analysis(
    [str(project_root / 'star_visibility.py')],
    pathex=[str(conda_bin)],
    binaries=binaries,
    datas=datas,
    hiddenimports=['_ssl', '_hashlib', 'ssl', 'astroquery', 'astroquery.simbad', 'astroquery.simbad.core', 'astroquery.query', 'astroquery.utils', 'astroquery.utils.commons', 'astroquery.utils.process_asyncs', 'astroquery.exceptions'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['PyQt6', 'PyQt5'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='StarVisibilityPlanner_simbad',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='StarVisibilityPlanner_simbad',
)
