from subprocess import Popen, PIPE



def run(cmd:str, cout:bool=False):
    process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = process.communicate()
    if stderr:
        print(stderr.decode("ascii"))
    if cout:
        print(stdout.decode("ascii"))