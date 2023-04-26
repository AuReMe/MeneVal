from meneval.meneval import *


def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--init', action='store_true', required=False, help='Init environment')
    parser.add_argument('--check', action='store_true', required=False, help='Check for required input files')
    parser.add_argument('--files', action='store_true', required=False, help='generate additional required input files')
    parser.add_argument('--blastp', action='store_true', required=False, help='Runs blastp step')
    parser.add_argument('--holobiont', action='store_true', required=False, help='Runs holobiont step')
    parser.add_argument('--aucome', action='store_true', required=False, help='Runs Aucome step')
    parser.add_argument('--group', type=str, required=False, metavar='group name', help='Group name for aucome step')
    parser.add_argument('--fill', action='store_true', required=False, help='Runs fill step')
    parser.add_argument('--workflow', action='store_true', required=False, help='Runs all steps')
    args = parser.parse_args()
    return args.init, args.check, args.files, args.blastp, args.holobiont, args.aucome, args.group, args.fill, \
        args.workflow


def main():
    INIT, CHECK, FILES_GENERATION, BLASTP_STEP, HOLOBIONT_STEP, AUCOME_STEP, GROUP, FILL_STEP, WORKFLOW = \
        get_command_line_args()

    if WORKFLOW:
        INIT = False
        CHECK = True
        FILES_GENERATION = True
        BLASTP_STEP = True
        HOLOBIONT_STEP = True
        AUCOME_STEP = True
        FILL_STEP = True

    # INITIALIZATION AND CHECK =========================================================================================
    if INIT:
        create_folders()
    if CHECK:
        check_required_files()

    # GENERATE MENECO FILES NEEDED =====================================================================================
    if FILES_GENERATION:
        generate_files()

    # RUN STEPS ========================================================================================================

    if BLASTP_STEP:
        run_step(1)

    if HOLOBIONT_STEP:
        run_step(2)

    if AUCOME_STEP:
        run_step(3, group=GROUP)

    if FILL_STEP:
        run_step(4)
        make_meneco_stats()


if __name__ == "__main__":
    main()
