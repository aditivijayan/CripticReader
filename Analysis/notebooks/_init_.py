import os, sys
script_dir = os.path.dirname(os.path.realpath('__file__'))
base_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
sys.path.append(base_dir)

sys.path.append(os.path.join(base_dir, 'config'))
sys.path.append(os.path.join(base_dir, 'lib'))
sys.path.append(os.path.join(base_dir, 'scripts'))                                             