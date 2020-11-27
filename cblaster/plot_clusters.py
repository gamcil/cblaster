
import tempfile
from pathlib import Path
import subprocess
import shutil
import logging


from cblaster.extract_clusters import extract_clusters


LOG = logging.getLogger(__name__)


def find_genbank_files(files):
    genbank_files = []
    for path in files:
        path_obj = Path(path)
        if path_obj.is_dir():
            genbank_files.extend(
                [str(po.resolve()) for po in path_obj.iterdir() if po.suffix in (".gbk", ".gb", ".genbank", ".gbff")])
        elif path_obj.suffix in (".gbk", ".gb", ".genbank", ".gbff"):
            genbank_files.append(str(path_obj.resolve()))
    return genbank_files


def run_clinker(genbank_files, allign_clusters, identity, plot_outfile, allignment_out):
    clinker_command = f"clinker {' '.join(genbank_files)} {'-na' if not allign_clusters else ''} -i {identity} " \
        f"-p {plot_outfile} {'-o' + allignment_out if allignment_out else ''}"
    process = subprocess.Popen(clinker_command, shell=True, stderr=subprocess.PIPE)

    # catch lines from stderr as they come in and raise when an error is encountered
    for line in process.stderr:
        print(line.decode(), end='')
        if "error" in line.decode().lower():
            process.kill()
            raise subprocess.CalledProcessError(-1, clinker_command)

    # make sure to kill the process after completion as per subprocess docs
    process.communicate()
    process.kill()


def plot_clusters(
    session=None,
    files=None,
    cluster_numbers=None,
    score_threshold=None,
    organisms=None,
    scaffolds=None,
    allign_clusters=False,
    identity=0.3,
    plot_outfile=None,
    allignment_out=None,
    cluster_out=None,
    prefix="",
):
    # if a session file is provided make genbank files first.
    remove_temp = False
    if session:
        if not cluster_out:
            cluster_out = tempfile.mkdtemp()
            remove_temp = True
        extract_clusters(
            session,
            cluster_out,
            file_format="genbank",
            prefix=prefix,
            cluster_numbers=cluster_numbers,
            score_threshold=score_threshold,
            organisms=organisms,
            scaffolds=scaffolds,
        )
        files = [cluster_out]

    # get only the genbank files present in directories and files
    genbank_files = find_genbank_files(files)

    try:
        run_clinker(genbank_files, allign_clusters, identity, plot_outfile, allignment_out)
    except subprocess.CalledProcessError:
        # make sure to remove the temp dir even when clinker crashes
        if remove_temp:
            shutil.rmtree(cluster_out)
        raise SystemExit
    # make sure to remove the temp dir
    if remove_temp:
        shutil.rmtree(cluster_out)
    LOG.info(f"Plot file can be found at {plot_outfile}")
    LOG.info("Done!")
