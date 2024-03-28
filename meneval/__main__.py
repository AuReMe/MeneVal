from meneval.meneval import *
import argparse


def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--init', action='store_true', required=False, help='Init environment')
    parser.add_argument('--check', action='store_true', required=False, help='Check for required input files')
    parser.add_argument('--files', action='store_true', required=False, help='generate additional required input files')
    parser.add_argument('--blastp', action='store_true', required=False, help='Runs blastp step')
    parser.add_argument('--enrich', type=str, required=False, metavar='group name', help='Group for enrichment step')
    parser.add_argument('--fill', action='store_true', required=False, help='Runs fill step')
    parser.add_argument('--workflow', action='store_true', required=False, help='Runs all steps')
    parser.add_argument('--exclude', action='store_true', required=False, help='Remove enrichment step reactions to network')
    args = parser.parse_args()
    return args.init, args.check, args.files, args.blastp, args.enrich, args.fill, args.workflow, args.exclude


def main():
    init, check, files_generation, blastp, enrich, fill, workflow, exclude_enrich = get_command_line_args()

    if workflow:
        check_required_files()
        generate_files()
        if check_step_required_files(BLASTP):
            run_step(BLASTP)
        group = GROUP_ALL
        if check_step_required_files(ENRICH, group):
            run_step(ENRICH, group)
        run_step(FILL)
        run_step(EXCLUDE_E)
        make_meneco_stats()

    # INITIALIZATION AND CHECK =========================================================================================
    if init:
        create_folders()
    if check:
        check_required_files()

    # GENERATE MENECO FILES NEEDED =====================================================================================
    if files_generation:
        generate_files()

    # RUN STEPS ========================================================================================================

    if blastp:
        if check_step_required_files(BLASTP):
            run_step(BLASTP)

    if enrich is not None:
        if enrich.upper() == GROUP_ALL:
            enrich = GROUP_ALL
        if check_step_required_files(ENRICH, enrich):
            run_step(ENRICH, group=enrich)

    if fill:
        run_step(FILL)
        make_meneco_stats()

    if exclude_enrich:
        run_step(EXCLUDE_E)


if __name__ == "__main__":
    main()
