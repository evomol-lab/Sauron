import sys
import re

with open('sauron.py', 'r') as f:
    content = f.read()

# Add rinpy import
if "import prody" in content:
    content = content.replace("import prody", "import prody\nexcept ImportError:\n    prody = None\n\ntry:\n    import rinpy\nexcept ImportError:\n    rinpy = None\n")
else:
    # Fallback if prody import is different
    content = content.replace("import argparse", "import argparse\ntry:\n    import rinpy\nexcept ImportError:\n    rinpy = None\n")


# Add rinpy logic in main
rinpy_logic = """
        if args.calc_method == 'rinpy':
            if rinpy is None:
                print("RinPy is not installed. Please install it using 'pip install rinpy'.")
                sys.exit(1)
            
            import shutil
            import os
            
            rinpy_in = f"{prefix}_rinpy_in"
            rinpy_out = f"{prefix}_rinpy_out"
            
            os.makedirs(rinpy_in, exist_ok=True)
            shutil.copy(input_file, os.path.join(rinpy_in, os.path.basename(input_file)))
            
            from rinpy import RINProcess
            proc = RINProcess(output_path=rinpy_out, input_path=rinpy_in, num_workers=1)
            try:
                proc.start_process()
            except Exception as e:
                print(f"RinPy calculation failed: {e}")
                sys.exit(1)
                
            print(f"RinPy calculation completed. Results are in {rinpy_out}")
            
            # Write Log
            import datetime
            with open(f"{prefix}_sauron.log", 'w') as logf:
                logf.write(f"Sauron RIN Calculation Log\\n")
                logf.write(f"Timestamp: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
                logf.write(f"Input File: {args.input}\\n")
                logf.write(f"Method: RINPY\\n")
                
            if os.path.exists(temp_pdb): os.remove(temp_pdb)
            if pqr_file != temp_pdb and os.path.exists(pqr_file): os.remove(pqr_file)
            shutil.rmtree(rinpy_in)
            return
"""

content = content.replace(
    "        if args.calc_method == 'insty':",
    rinpy_logic + "\n        if args.calc_method == 'insty':"
)

with open('sauron.py', 'w') as f:
    f.write(content)

