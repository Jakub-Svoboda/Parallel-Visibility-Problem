from subprocess import Popen, PIPE
from argparse import ArgumentParser
import math
import random
import itertools

FAIL = '\033[91m\033[1m'
GREEN = '\033[32m\033[1m'
ENDC = '\033[0m'
BOLD = '\033[1m'


def main():
    parser = ArgumentParser()
    parser.add_argument("-m", "--mpi", type=str)
    parser.add_argument("-e", "--executable", type=str, default="vid")
    parser.add_argument("-l", "--limit", type=int, default=15)

    argv = parser.parse_args()
    if not argv.executable:
        return

    altitudes = [random.randint(1, 1024)]

    args = [argv.mpi, '-np', '', argv.executable, '']
    for inputSize in range(1, 31):
        print(f'[Input size] {inputSize}')
        # permutations = itertools.permutations(altitudes)

        for i in range(0, min(argv.limit, math.factorial(inputSize))):
            if i == 0:
                altitudes.sort()
            elif i == 1:
                altitudes.sort(reverse=True)
            else:
                random.shuffle(altitudes)

            args[4] = ','.join(map(str, altitudes))
            print(f'[Input string] {args[4]}')
            angles = [-math.inf] + [math.atan((x - altitudes[0]) / float(i))
                    for i, x in enumerate(altitudes[1:], 1)]

            max_prescan = [-math.inf] + list(itertools.accumulate(angles, max))
            max_prescan.pop()
            visibilities = ["v" if (angle > max_prev_angle) else "u"
                            for (angle, max_prev_angle) in zip(angles, max_prescan)]
            visibilities[0] = "_"

            ref_output = ",".join(visibilities)
            #for procCount in range(1, inputSize + 1):
            #args[2] = str(procCount)
            
            myargs = ["./nocomptest.sh", args[4]]
            print(myargs)
            process = Popen(myargs, stdout=PIPE, stderr=PIPE)
            output = process.communicate()[0].decode('utf-8').rstrip('\n')
            if output != ref_output:
                print(f"{FAIL}FAIL: {ENDC}{FAIL}{output}{ENDC} != {GREEN}{ref_output}{ENDC}")
                print(f'[Processor count] {procCount}')

        altitudes.append(random.randint(1, 1024))


if __name__ == "__main__":
    main()
