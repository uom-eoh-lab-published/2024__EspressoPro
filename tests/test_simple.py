#!/usr/bin/env python3
import sys
from pathlib import Path

# Add src to path
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

try:
    import core
    
    get_default_models_path = core.get_default_models_path
    get_default_data_path = core.get_default_data_path
    
    print('✅ Functions imported successfully')
    
    # Test models path
    try:
        models_path = get_default_models_path()
        print(f'🤖 Models path: {models_path}')
        print(f'   Exists: {models_path.exists()}')
    except Exception as e:
        print(f'❌ Models path error: {e}')
    
    # Test data path
    try:
        data_path = get_default_data_path()
        print(f'📊 Data path: {data_path}')
        print(f'   Exists: {data_path.exists()}')
    except Exception as e:
        print(f'❌ Data path error: {e}')
        
except ImportError as e:
    print(f'Import error: {e}')
