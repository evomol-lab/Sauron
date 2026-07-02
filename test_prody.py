import sys
sys.path.insert(0, '/home/jpmslima/coding/Sauron')
try:
    import prody
    print("Top level prody OK")
except Exception as e:
    print(f"Top level prody error: {e}")

def test_local_import():
    try:
        import prody
        print("Local prody OK")
    except ImportError as e:
        print(f"Local prody error: {e}")

test_local_import()
