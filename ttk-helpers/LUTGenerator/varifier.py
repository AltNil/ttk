import sys
import subprocess

import config

TOO_BIG_SIZE = 150
def getCaseId(fromid, size):
    factor = 3
    dimension = [size, size, size]
    coords = []
    intermediate = fromid
    dimProd = 1
    for i in dimension[:-1]:
        dimProd = dimProd * int(i)
    for i in dimension[::-1]:
        coord = int(intermediate/dimProd)
        intermediate = intermediate -coord*dimProd
        dimProd = dimProd/i
        coords.append(coord)
    coords.reverse()
    res = 0
    for i in range(len(dimension)):
        if int(coords[i])==0:
            res = res+factor**i
        if int(coords[i])==dimension[i]-1:
            res = res+2*(factor**i)
    return res


def main():
    tooBig = False
    if len(sys.argv) < 2:
        print("Please provide a input file as argument")
        return
    output_format = 0
    if len(sys.argv) >= 3:
        if sys.argv[2] == "Short":
            output_format = 1
        elif sys.argv[2] == "Long":
            output_format = 0
        elif sys.argv[2] == "ExportShort":
            output_format = 3
        elif sys.argv[2] == "ExportLong":
            output_format = 2
        elif sys.argv[2] == "ExtremKurz":
            output_format = -1
    else:
        print("Output format options not used, available are:\n'Short' 'Long' 'ExportShort' 'ExportLong'")
    if output_format >= 0:
        print("outputformat set to ", output_format)
        print("Starting validation...")
    res = None
    testCase = sys.argv[1]
    testSize = int(sys.argv[1][-8:-5])+1
    try:
        print("Executing testcase " + testCase + " with size " + str(testSize))
        res = subprocess.check_output(
            [config.ttk_path, "--exit", testCase]).decode(
            sys.stdout.encoding)
    except:
        print("Something somewhere went wrong.")
    if output_format >= 0:
        print("Parsing the output...")
        tooBig = True
    elif testSize > TOO_BIG_SIZE:
        print("Inputsize of the grid is too big, only handle timing...")
        tooBig = True
    batch = 0
    case1 = {}
    c1Timing = "Dummytiming fromTmie[0000s]"
    case2 = {}
    c2Timing = "Dummytiming fromTmie[0000s]"
    seperator = "========"
    try:
        for line in res.splitlines():
            if batch == 0:
                if seperator in line:  # Begin of first test
                    batch = 1
                continue
            elif "Completed" in line:  # End of tests
                if batch == 1:
                    c1Timing = line
                    batch = 2
                    continue
                elif batch == 2:
                    c2Timing = line
                    break

            if seperator in line or tooBig:
                continue
            lineSplit = line.split(":")
            if batch == 1:
                case1[lineSplit[0]] = lineSplit[1] if len(lineSplit) == 2 else ""
            elif batch == 2:
                case2[lineSplit[0]] = lineSplit[1] if len(lineSplit) == 2 else ""
        # Start comparison
        if output_format >= 0:
            print("Checking the Results...")
        correctness = case1 == case2
    except Exception as e:
        print("Something unexpected happened")
        print(e)
        correctness = True
    if output_format >= 0:
        print("Result of comparison:", correctness if not tooBig else " Not possible")
    if not correctness:
        count = 0
        print("Here are the different lines")
        for key in case1.keys():
            if case1[key] == case2[key]:
                continue
            prepKey = key[key.find(" ")+5:]
            idStart = prepKey.find("[")
            fromId = (prepKey[idStart+1:prepKey.find("]")])
            caseId = -1
            prepKey = prepKey[0:idStart]
            try:
                caseId = getCaseId(int(fromId), testSize)
            except:
                print("could not transform to int")
            sep = ","
            c1Split = case1[key].split(sep)
            c2Split = case2[key].split(sep)
            delta = ["0"] * len(c1Split)
            diff = []
            for i in range(len(c1Split)):
                if c1Split[i] != c2Split[i]:
                    diff.append(i)
            for i in diff:
                highlight = "\x1b[0;31;40m"
                lowlight = "\x1b[0m"
                if output_format == 2 or output_format == 3:
                    highlight = ""
                    lowlight = ""
                try:
                    v1 = c1Split[i]
                    v2 = c2Split[i]
                    if i == 0:
                        v1 = v1[v1.find("[")+1:]
                        v2 = v2[v2.find("[")+1:]
                    delta[i] =highlight +  str(int(v1)-int(v2)) + lowlight
                except:
                    print("Geht halt einfach nicht...")
                if i == 0:
                    c1Split[i] = highlight + c1Split[i][c1Split[0].find("[")+1:] + lowlight
                    c2Split[i] = highlight + c2Split[i][c2Split[0].find("[")+1:] + lowlight
                else:
                    c1Split[i] = highlight + c1Split[i] + lowlight
                    c2Split[i] = highlight + c2Split[i] + lowlight
            case1[key] = "[" + sep.join(c1Split)
            case2[key] = "[" + sep.join(c2Split)
            deltaString = "[" + sep.join(delta) + "]"
            if output_format%2 == 0:
                print(prepKey, "(id =",  fromId, "with caseId =", str(caseId) + ") diffs at ", diff, "\n", case1[key], "\n", case2[key], "\n", deltaString)
            elif output_format%2 == 1:
                print(prepKey, "(id =", fromId, "with caseId =", str(caseId) + ") diffs at ", diff)
            count = count +1
        print(count, " values were different")
    print("Timings:\n", c1Timing[c1Timing.find("Completed"):], "\n", c2Timing[c1Timing.find("Completed"):])
    time1 = float(
        c1Timing[c1Timing.find("[", c1Timing.find("[", c1Timing.find("Completed")) + 3) + 1:c1Timing.find("|") - 1])
    time2 = float(
        c2Timing[c2Timing.find("[", c2Timing.find("[", c2Timing.find("Completed")) + 3) + 1:c2Timing.find("|") - 1])
    print("Speedup:", time1 / time2 if time2 != 0 else "infinity")


main()
