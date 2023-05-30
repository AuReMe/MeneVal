from meneval.meneval import *
import argparse


def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--init', action='store_true', required=False, help='Init environment')
    parser.add_argument('--check', action='store_true', required=False, help='Check for required input files')
    parser.add_argument('--files', action='store_true', required=False, help='generate additional required input files')
    parser.add_argument('--blastp', action='store_true', required=False, help='Runs blastp step')
    parser.add_argument('--enrich', type=str, required=False, metavar='group name', help='Group name for aucome step')
    parser.add_argument('--fill', action='store_true', required=False, help='Runs fill step')
    parser.add_argument('--workflow', action='store_true', required=False, help='Runs all steps')
    args = parser.parse_args()
    return args.init, args.check, args.files, args.blastp, args.enrich, args.fill, args.workflow


def main():
    init, check, files_generation, blastp, enrich, fill, workflow = get_command_line_args()

    if workflow:
        init = False
        check = True
        files_generation = True
        blastp = True
        fill = True

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
        if check_step_required_files(ENRICH, enrich):
            run_step(ENRICH, group=enrich)

    if fill:
        run_step(FILL)
        make_meneco_stats()


if __name__ == "__main__":
    main()
