with open('sauron.py', 'r') as f:
    content = f.read()

content = content.replace(
    'parser.add_argument("--chains", type=str, help="Comma-separated list of chains to calculate (e.g. A,B,C)")',
    'parser.add_argument("--chains", type=str, help="Comma-separated list of chains to calculate (e.g. A,B,C)")\n    parser.add_argument("--energy", action="store_true", help="Calculate Interaction Energy (PENs). Requires a PQR file or --add-h.")\n    parser.add_argument("--vdw-energy", action="store_true", help="Include generic Lennard-Jones VDW term in energy calculation. If not set, only Electrostatics is computed.")'
)

content = content.replace(
    'edges_df, nodes = calculate_rin(struct_to_calc, strict_angle=args.strict_angle, remove_multiples=args.remove_multiples, model_num=model_num, calc_method=args.calc_method, pqr_file=pqr_file)',
    'edges_df, nodes = calculate_rin(struct_to_calc, strict_angle=args.strict_angle, remove_multiples=args.remove_multiples, model_num=model_num, calc_method=args.calc_method, pqr_file=pqr_file, calc_energy=args.energy, include_vdw=args.vdw_energy)'
)

with open('sauron.py', 'w') as f:
    f.write(content)

