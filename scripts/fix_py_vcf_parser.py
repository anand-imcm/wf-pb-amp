import sys, os
import shutil

"""
To address an issue in the default parser provided by the PyVCF package,
causing the following error message: 'TypeError: "quotechar" must be a 1-character string.'
This error disrupts the functionality of the VCFCons.py script in producing VCF output files.

Solution:
The problem arises from the default value of the 'quotechar' parameter, which cannot be set to None
within the PyVCF parser. Instead, it should be set to '"' (double quotation mark),
which aligns with the expected default value of the csv.writer() function.
"""

script_dir = os.path.dirname(os.path.abspath(__file__))
new_parser_path = os.path.join(script_dir,"parser.py")

try:
    import vcf
    for path in sys.path:
        parser_path = os.path.join(path, "vcf","parser.py")
        if "site-packages" in path and os.path.exists(parser_path):
            shutil.copy(new_parser_path, parser_path)
            print (f" ---> Fixed : {parser_path}")
except ModuleNotFoundError:
    print(f" ---> Warning : pyvcf module can't be used in this pipeline")

