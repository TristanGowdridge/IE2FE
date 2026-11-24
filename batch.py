import os
import re
import subprocess

import click
from tqdm import tqdm


def validate_virual_path(ctx, param, value):
    if value != None:
        suffix_re = re.compile(r"^[a-zA-Z0-9_-]+(\.[a-zA-Z]+)?$")
        if not suffix_re.fullmatch(value):
            raise click.BadParameter(
                "format must be alpha numeric with an optional _ or -"
            )
    return value


virtual_path = dict(
    required=True, type=click.UNPROCESSED, callback=validate_virual_path
)


@click.command()
@click.option(
    "-d",
    "--base-directory",
    required=True,
    type=click.Path(
        exists=True,
        readable=True,
        writable=True,
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
    ),
    multiple=True,
)
@click.option("-ie", "--ie-directory", **virtual_path)
@click.option("-fe", "--fe-directory", **virtual_path)
@click.option("-l", "--loadcase-directory", **virtual_path)
@click.option("-r", "--result-directory", **virtual_path)
@click.option("-f", "--feature-directory", **virtual_path)
@click.option("-p", "--parameters-file", **virtual_path)
@click.option("-log", "--log-file", default=None, **virtual_path)
@click.option("--verbose", required=False, is_flag=True)
@click.option("--limit", required=False, type=int, default=0)
@click.option("--dispatch-clear", required=False, type=int, default=10)
def batch(
    base_directory: str,
    ie_directory: str,
    fe_directory: str,
    loadcase_directory: str,
    result_directory: str,
    feature_directory: str,
    parameters_file: str,
    log_file: str,
    verbose: bool,
    limit: int,
    dispatch_clear: int,
) -> None:

    # create run parameters
    runs = []
    for directory in base_directory:
        parameters = {
            f"-{option}": os.path.join(directory, value)
            for option, value in [
                ("ie", ie_directory),
                ("fe", fe_directory),
                ("l", loadcase_directory),
                ("r", result_directory),
                ("f", feature_directory),
                ("p", parameters_file),
                ("log", log_file),
            ]
            if value is not None
        } | {
            f"--{option}": str(value)
            for option, value in [
                ("limit", limit),
                ("dispatch-clear", dispatch_clear),
            ]
        }

        runs.append(parameters)

        print(f"Created the following parameters for {directory}")
        for key, value in parameters.items():
            print(f"\t{key}: {value}")

    # run the batches
    for options in tqdm(runs, desc="Runs"):
        run_options = [item for option in options.items() for item in option]
        run = subprocess.run(
            (["python", "run.py"] + run_options + (["--verbose"] if verbose else []))
        )
        if run.returncode != 0:
            print(
                f"Error running batch, return code: {run.returncode} "
                f"with options: {run_options}"
            )


if __name__ == "__main__":
    batch()
