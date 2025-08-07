#!/usr/bin/env python3
"""
Test automatic model availability checks in annotation functions
"""

def test_model_availability_checks():
    """Test that annotation functions include model availability checks"""
    
    # Test imports work
    try:
        import sys
        sys.path.insert(0, 'src')
        from annotation import annotate_anndata, annotate_counts_matrix
        from missionbio import annotate_missionbio_sample
        from core import ensure_models_available, download_models
        print("✅ All annotation functions imported successfully")
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    
    # Test that functions have the ensure_models_available check
    import inspect
    
    # Check annotate_anndata
    source = inspect.getsource(annotate_anndata)
    if "ensure_models_available" in source and "Checking model availability" in source:
        print("✅ annotate_anndata includes model availability check")
    else:
        print("❌ annotate_anndata missing model availability check")
        return False
    
    # Check annotate_counts_matrix  
    source = inspect.getsource(annotate_counts_matrix)
    if "ensure_models_available" in source and "Checking model availability" in source:
        print("✅ annotate_counts_matrix includes model availability check")
    else:
        print("❌ annotate_counts_matrix missing model availability check")
        return False
    
    # Check annotate_missionbio_sample
    source = inspect.getsource(annotate_missionbio_sample)
    if "ensure_models_available" in source and "Checking model availability" in source:
        print("✅ annotate_missionbio_sample includes model availability check")
    else:
        print("❌ annotate_missionbio_sample missing model availability check")
        return False
    
    print("✅ All annotation functions include automatic model availability checks!")
    return True

if __name__ == "__main__":
    print("Testing EspressoPro model availability checks...")
    print("=" * 60)
    
    success = test_model_availability_checks()
    
    if success:
        print("\n🎉 All tests passed! Automatic model downloading is ready.")
    else:
        print("\n❌ Some tests failed. Please check the implementation.")
