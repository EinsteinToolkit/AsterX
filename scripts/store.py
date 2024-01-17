'''
This file stores logs for future use in the records folder
'''
import glob
import os
import shutil
import subprocess
import sys
import configparser

def copy_tests(test_dir,version,procs):
    '''
        This function copies individual test logs and diffs.
        It takes the directory of the logs, the version number 
        and then number of proceses.
    '''
    dst=f"./records/version_{version}/sim_{version}_{procs}"
    dirlist=os.listdir(test_dir)
    if not os.path.isdir(dst):
        os.mkdir(dst)
    for dir in dirlist:
        if(not os.path.isdir(dst+"/"+dir)):
            os.mkdir(dst+"/"+dir)
        diffs=[x.split("/")[-1] for x in glob.glob(test_dir+"/"+dir+"/*.diffs")]
        logs=[x.split("/")[-1] for x in glob.glob(test_dir+"/"+dir+"/*.log")]
        for log in logs:
            shutil.copy(test_dir+"/"+dir+"/"+log,dst+"/"+dir+"/"+log)
        for diff in diffs:
            shutil.copy(test_dir+"/"+dir+"/"+diff,dst+"/"+dir+"/"+diff)

    #shutil.copytree(test_dir,dst,)


def copy_logs(test_dir,version):
    '''
        This copies the test logs for future use
    '''
    dst=f"./records/version_{version}/"
    log=f"{test_dir}/summary.log"
    # TODO: move into separate function and return Python object
    config=configparser.ConfigParser()
    config.read(f"{test_dir}/../../SIMFACTORY/properties.ini")
    procs=config['properties']['procs']
    numthreads=config['properties']['numthreads']
    build=f"build__{procs}_{numthreads}_{version}.log"
    shutil.copy(log,dst+build)


def copy_index(version):
    '''
        This copies the old html files showing test
        results for future use
    '''
    dst=f"./docs/index_{version}.html"
    index="./docs/index.html"
    if os.path.exists(index):
        shutil.copy(index,dst)

def copy_compile_log(version):
    '''
        This copies the compilation logs for future use
    '''
    GITHUB_WORKSPACE=os.environ["GITHUB_WORKSPACE"]
    dst=f"./records/version_{version}/build_{version}.log"
    # This must match the build command in `build.sh`
    build=f"{GITHUB_WORKSPACE}/../workspace/Cactus/build.log"
    shutil.copy(build,dst)

def store_commit_id(version):
    '''
        This stores the current git HEAD hash for future use
    '''
    global REPO
    # Keep trailing newline with the commit id
    id=subprocess.check_output(["git", "--git-dir", f"{REPO}/.git", "rev-parse", "HEAD"])
    dst=f"./records/version_{version}/id.txt"
    with open(dst,'wb') as f:
        f.write(id)

def get_version():
    '''
        This checks the version of the current build
        by looking at the file names from old builds.
    '''
    current=0
    build_records=glob.glob("./records/version_*")
    builds=[int(x.split("_")[-1].split(".")[0]) for x in build_records]
    try:
        current_build=max(builds)
    except ValueError:
        current_build = 0
    return current_build

def store_version(next_build):
    '''
        This stores the version of the current build
        in the list of build numbers file.
    '''
    with open("./docs/version.txt",'a') as vers:
        vers.write(f"{next_build}\n")

def get_commit_id(version):
    '''
        Returns the code commit id that this version corresponds to.
    '''
    try:
        with open(f"./records/version_{version}/id.txt", "r") as fh:
            id = fh.readline().strip()
    except FileNotFoundError:
        id = "0"
    return id

if __name__ == "__main__":
    # FIXME: this is quite bad, use some better argparse
    REPO = sys.argv[1]
    dir1 = sys.argv[2]
    dir2 = sys.argv[3]
    version=get_version()+1
    store_version(version)
    os.mkdir(f"./records/version_{version}/")
    copy_compile_log(version)
    copy_logs(dir1,version)
    copy_logs(dir2,version)
    copy_tests(dir1,version,1)
    copy_tests(dir2,version,2)
    store_commit_id(version)
