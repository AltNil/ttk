import sys
import subprocess
import re

import config

factor = 3
dims = 3
dimension = []


neighbor_list = []
#
# Calculates the coordinates of a point by its id
#
def calculatePCoords(pid):
    coords = []
    intermediate = pid
    dimProd = 1
    for i in dimension[:-1]:
        dimProd = dimProd * int(i)
    for i in dimension[::-1]:
        coord = int(intermediate/dimProd)
        intermediate = intermediate -coord*dimProd
        dimProd = dimProd/i
        coords.append(coord)
    coords.reverse()
    return coords


#
# calculates the caseId given a point id
#
def caseCalculation(pid):
    coords = calculatePCoords(pid)
    res = 0
    for i in range(len(dimension)):
        if int(coords[i])==0:
            res = res+factor**i
        if int(coords[i])==dimension[i]-1:
            res = res+2*(factor**i)
    return res


#
# Taking as input an list of arrays in which the position 1 referes to form where or
# the nNeighbors and the position 2 contains the information where or how many
#
def prepareInputTag(listInput, nTag, tag):
    nTagOut = []
    tagOut = [None]*(len(listInput) - 1)
    for line in listInput:
        op = line[2].find("[")
        en = line[2].find("]")
        idSplit = re.split(",",line[2][op+1:en-1])
        if nTag in line[1]:
            nTagOut = idSplit
        elif tag in line[1]:
            op2 = line[1].find("s[")
            en2 = line[1].find("]")
            fromId = int(line[1][op2+2:en2])
            tagOut[fromId]=idSplit
        else:
            continue
    return nTagOut, tagOut
            

def calculateNTable(nList):
    lut = [0]*(factor**dims)
    for i in range(len(nList)):
        case = caseCalculation(i)
        lut[case] = nList[i]
    return lut


def intArrayString(arrayName ,nList):
    out = "const int "+arrayName+"[" + str(len(nList)) + "] = {"
    for i in nList:
        out = out + str(i) + ","
    out = out[0:-1] + "};"
    return out


def handleNeighbors(inputList):
    resString = ""
    nNeighbor, neighbors = prepareInputTag(inputList, "nNeighbors", "neighbor")
    # prepare outputVariables
    lutNeighbors = [None] * (factor ** dims)
    lutnNeighbors = calculateNTable(nNeighbor)
    lutEdgeDirections = [None] * (factor ** dims)

    #
    # calculating the coordinate offsets
    #
    for i in range(len(neighbors)):
        toIDs = neighbors[i]
        if (toIDs is None):
            continue
        # toIDs = line.split(",")
        fromCoords = calculatePCoords(i)
        caseID = caseCalculation(i)

        deltas = [[0, 0, 0]] * len(toIDs)
        directions = [0] * len(toIDs)
        j = 0
        for to in toIDs:
            num = int(to)
            newCoords = calculatePCoords(num)
            offset = [0, 0, 0]
            for i in range(3):
                offset[i] = newCoords[i] - fromCoords[i]
            deltas[j] = offset
            directions[j] = offset[0]+1+3*(offset[1]+1)+9*(offset[2]+1)
            j = j + 1
        lutNeighbors[caseID] = deltas
        lutEdgeDirections[caseID] = directions
    resString = resString + "\n\n" +  intArrayString("lutNumNeighbor3d",lutnNeighbors)

    # formatting the output
    outNeighborOffset = "const std::vector<std::vector<std::array<int,3>>> lutNeighborOffset = {"
    for li in lutNeighbors:
        if li is None:
            li = "{}"
        listr = str(li).replace("[", "{").replace("]", "}")
        outNeighborOffset = outNeighborOffset + "\n\t" + listr + ","
    outNeighborOffset = outNeighborOffset[0:-1] + "\n};"
    resString = resString + "\n\n" + outNeighborOffset

    outEdgeDirections = "const std::vector<std::vector<int>> lutEdgeDirection3d = {"
    for li in lutEdgeDirections:
        if li is None:
            li = "{}"
        listr = str(li).replace("[","{").replace("]","}")
        outEdgeDirections = outEdgeDirections + "\n\t" + listr + ","
    outEdgeDirections = outEdgeDirections[0:-1] + "\n};"
    resString = resString + "\n\n" + outEdgeDirections

    return resString


def threeSolver(solution, offset, x, y, z):
    res = [0, 0, 0]
    calculation = solution - offset
    while z < calculation:
        res[2] += 1
        calculation -= z
    while y < calculation:
        res[1] += 1
        calculation -= y
    while x <= calculation:
        res[0] += 1
        calculation -= x
    if calculation != 0:
        return [calculation]
    if calculation == 0:
        return res
    return "was soll das?"


def prepareNeighborList(tetra = True):
    global neighbor_list
    if tetra:
        neighbor_list = [[-1, 0, 0],[0, -1, -1],[0, 0, 0],[0, -1, 0],[-1, -1, 0],[0, 0, -1],[-1, 0, -1],[-1, -1, -1]]
    else:
        neighbor_list = [[0, 0, 0],[-1, 0, 0],[1, -1, -1],[0, -1, -1],[0, -1, 0],[1, -1, 0],[0, 0, -1],[1, 0, -1]]
    print("Preperation complete\n", neighbor_list)


def optionSolver(solution, offset, x, y, z, cx, cy, cz):
    options = neighbor_list
    for option in options:
        op = [option[0]+cx, option[1]+cy, option[2]+cz]
        possible_solution = op[0]*x + op[1]*y + op[2]*z + offset
        if possible_solution == solution:
            return op
    return [-1]


def getDirectionAndOffset(fromId, triangleId:int):
    global dimension
    dx = dimension[0]
    dy = dimension[1]
    dz = dimension[2]

    plane1 = (dx - 1) * 2 * (dy - 1) * dz
    plane2 = plane1 + (dx - 1) * 2 * (dy) * (dz - 1)
    plane3 = plane2 + (dx * (dy - 1) * (dz - 1) * 2)
    plane4 = plane3 + (dx - 1) * (dy - 1) * (dz - 1) * 2
    plane5 = plane4 + (dx - 1) * (dy - 1) * (dz-1) * 2

    dxy = dx*dy
    cz = int(fromId/dxy)
    xy = fromId-cz*dxy
    cy = int(xy/dx)
    cx = xy-cy*dx
    current = [cx, cy, cz]

    if int(0) <= triangleId < plane1: # 0 -> 23
        x = 2
        y = 2*(dx-1)
        z = 2*(dx-1)*(dy-1)
        offset = 0
        possibility1 = optionSolver(triangleId, offset, x, y, z, cx, cy, cz)
        possibility2 = optionSolver(triangleId, offset-1, x, y, z, cx, cy, cz)
        if len(possibility2) == 1:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility1[i]-current[i]
            return out, 446
        else:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility2[i]-current[i]
            return out, 447
    elif plane1 <= triangleId < plane2: # 24 -> 47
        x = 2
        y = 2*(dx-1)
        z = 2*(dx-1)*(dy)
        offset = plane1
        possibility1 = optionSolver(triangleId, offset, x, y, z, cx, cy, cz)
        possibility2 = optionSolver(triangleId, offset-1, x, y, z, cx, cy, cz)
        if len(possibility2) == 1:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility1[i]-current[i]
            return out, 608
        else:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility2[i]-current[i]
            return out, 609
    elif plane2 <= triangleId < plane3: # 48 -> 71
        x = 2
        y = 2*dx
        z = 2*(dx)*(dy-1)
        offset = plane2
        possibility1 = optionSolver(triangleId, offset, x, y, z, cx, cy, cz)
        possibility2 = optionSolver(triangleId, offset+1, x, y, z, cx, cy, cz)
        if len(possibility2) == 1:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility1[i]-current[i]
            return out, 691
        else:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility2[i]-current[i]
            return out, 697
    elif plane3 <= triangleId < plane4: # 72 -> 87
        x = 2
        y = 2*(dx-1)
        z = 2*(dx-1)*(dy-1)
        offset = plane3-2
        possibility1 = optionSolver(triangleId, offset, x, y, z, cx, cy, cz)
        possibility2 = optionSolver(triangleId, offset+1, x, y, z, cx, cy, cz)
        if len(possibility2) == 1:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility1[i]-current[i]
            return out, 659
        else:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility2[i]-current[i]
            return out, 669
    elif plane4 <= triangleId < plane5: # 88 -> 103
        x = 2
        y = 2*(dx-1)
        z = 2*(dx-1)*(dy-1)
        offset = plane4
        possibility1 = optionSolver(triangleId, offset, x, y, z, cx, cy, cz)
        possibility2 = optionSolver(triangleId, offset-1, x, y, z, cx, cy, cz)
        if len(possibility2) == 1:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility1[i]-current[i]
            return out, 689
        else:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility2[i]-current[i]
            return out, 699
    elif plane5 <= triangleId: # 104 ->
        x = 2
        y = 2*(dx-1)
        z = 2*(dx-1)*(dy-1)
        offset = plane5-2
        possibility1 = optionSolver(triangleId, offset, x, y, z, cx, cy, cz)
        possibility2 = optionSolver(triangleId, offset+1, x, y, z, cx, cy, cz)
        if len(possibility2) == 1:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility1[i]-current[i]
            return out, 717
        else:
            out = [0,0,0]
            for i in range(3):
                out[i]=possibility2[i]-current[i]
            return out, 724
    return 0, 0


def handleTriangles(listInput):
    prepareNeighborList(False)
    outString = "\n\n"
    nTriangles, triangle = prepareInputTag(listInput,"nTriangles", "triangle")
    lutNTriangles = calculateNTable(nTriangles)
    outString += intArrayString("lutVertexTriangles3d",lutNTriangles)

    #
    # Calculate the neighbors by offset and directionId
    #
    
    # Case distinction by number range
    lutTriOffsets = [None] * (factor ** dims)
    lutTriDirections = [0] * (factor ** dims)
    i = 0
    for arr in triangle:
        deltas = [[0, 0, 0]] * len(arr)
        directions = [0] * len(arr)

        for j in range(len(arr)):
            id = arr[j]
            offset, direction = getDirectionAndOffset(i,int(id))
            deltas[j] = offset
            directions[j] = direction
        caseId = caseCalculation(i)
        lutTriDirections[caseId] = directions
        lutTriOffsets[caseId] = deltas
        i +=1
    # resString = intArrayString("lutTriangleDirections3d", lutTriDirections)
    outTriangleDirections = ("const std::vector<std::vector<int>> lutTriangleDirection3d = {")
    for li in lutTriDirections:
        if li is None:
            li = "{}"
        liStr = str(li).replace("[","{").replace("]","}")
        outTriangleDirections = outTriangleDirections + "\n\t" + liStr + ","
    outTriangleDirections = outTriangleDirections[0:-1] + "\n};"
    outString = outString + "\n\n" + outTriangleDirections

    outTriangleOffset = ("const std::vector<std::vector<std::array<int,3>>> lutTriangleOffset3d = {")
    for li in lutTriOffsets:
        if li is None:
            li = "{}"
        listr = str(li).replace("[", "{").replace("]", "}")
        outTriangleOffset = outTriangleOffset + "\n\t" + listr + ","
    outTriangleOffset = outTriangleOffset[0:-1] + "\n};"
    outString = outString + "\n\n" + outTriangleOffset
    return outString


def handleTetrahedrons(listInput):
    prepareNeighborList(True)
    outString = "\n\n"
    nTetrahedrons, tetrahedrons = prepareInputTag(listInput,"nStars", "star")
    lutNTriangles = calculateNTable(nTetrahedrons)
    outString += intArrayString("lutVertexTetrahedrons3d",lutNTriangles)

    dx = dimension[0]
    dy = dimension[1]
    dz = dimension[2]

    ### Calculating tetrahedron directions
    lutTetraDirections = [0] * (factor ** dims)
    i = 0
    lutTetraOffsets = [None] * (factor ** dims)
    for arr in tetrahedrons:
        directions = [0] * len(arr)
        offset = []

        dxy = dx * dy
        cz = int(i / dxy)
        xy = i - cz * dxy
        cy = int(xy / dx)
        cx = xy - cy * dx

        plane1 = dx-1 * 2 * dy-1 * dz
        plane2 = plane1 + dx-1 * 2 * (dy) * dz-1
        plane3 = plane2 + (dx * dy-1 * dz-1 * 2)
        plane4 = plane3 + dx-1 * dy-1 * dz-1 * 2
        plane5 = plane4 + dx-1 * dy-1 * dz-1 * 2

        plane = [0,plane1,plane2,plane3,plane4,plane5]

        for j in range(len(arr)):
            id = arr[j]
            solverSolution = threeSolver(int(id), 0, 6,6*(dx-1),6*(dx-1)*(dy-1))
            direction = 0 if len(solverSolution) == 3 else solverSolution[0]
            directions[j] = direction
            offset_option = []
            for off in plane:
                total_offset = off + direction
                offset_option = optionSolver(int(id), total_offset, 6,6*(dx-1),6*(dx-1)*(dy-1), cx, cy, cz)
                if len(offset_option) == 3:
                    offset_option[0] = offset_option[0]-cx
                    offset_option[1] = offset_option[1]-cy
                    offset_option[2] = offset_option[2]-cz
                    break
            offset.append(offset_option)
        caseId = caseCalculation(i)
        lutTetraOffsets[caseId] = offset
        lutTetraDirections[caseId] = directions
        i +=1


    outTetraDirections = ("const std::vector<std::vector<int>> lutTetrahedronDirection3d = {")
    for li in lutTetraDirections:
        if li is None:
            li = "{}"
        liStr = str(li).replace("[","{").replace("]","}")
        outTetraDirections = outTetraDirections + "\n\t" + liStr + ","
    outTetraDirections = outTetraDirections[0:-1] + "\n};"
    outString = outString + "\n\n" + outTetraDirections

    outTriangleOffset = ("const std::vector<std::vector<std::array<int,3>>> lutTetraOffset3d = {")
    for li in lutTetraOffsets:
        if li is None:
            li = "{}"
        listr = str(li).replace("[", "{").replace("]", "}")
        outTriangleOffset = outTriangleOffset + "\n\t" + listr + ","
    outTriangleOffset = outTriangleOffset[0:-1] + "\n};"
    outString = outString + "\n\n" + outTriangleOffset
    return outString


def processList(l, tag):
    lookupTag = tag[1:-1]
    resString = ""
    if "neighbor" in lookupTag:
        resString += handleNeighbors(l)
    if "triangle" in lookupTag:
        resString += handleTriangles(l)
    if "star" in lookupTag:
        resString += handleTetrahedrons(l)
    return resString


def getNumLists(tag):
    lookupTag = tag[1:-1]
    if (lookupTag == "neighbor"):
        return 2
    return 0


def run():  
    if len(sys.argv)!=4:
        print("Needs 3 parameters:\n1. Path to test csv file containing testfile;gridSize\n2. #cases in each dimension\n3. #dimensions")
        return
    global factor
    global dims
    global dimension    
    tests = sys.argv[1]
    factor = int(sys.argv[2])
    dims = int(sys.argv[3])
    
    # read the CSV
    with open(tests) as f:
        lines = f.readlines()

    
    # execute tests
    for line in lines:
        test = line.split(";")
        if len(test) == 3:
            continue
        print("\t processing ", test[0])
        dimension = list(map(int,test[1].split(",")))
        try:
            res = subprocess.check_output([config.ttk_path, "--exit", test[0]]).decode(sys.stdout.encoding)
        except:
            print("an error occurred")
        started = False
        relevant = []
        currentTag = ""
        tagList = []
        tagSizes = []
        resultCode = ""

        # Filtering teh output for relevant lines
        for line in res.splitlines():
            if "nVertices" in str(line):
                started = True
            elif started and "========================================" in str(line):
                resultCode += processList(relevant, currentTag)
                break
            elif started:
                interestingStart = line.find(" ")
                interestingEnd = line.find("[", interestingStart+5)
                optionalEnd = line.find(":", interestingStart+5)
                tag = line[interestingStart+5:interestingEnd-1 if interestingEnd < optionalEnd else optionalEnd]
                if tag not in currentTag:
                    resultCode += processList(relevant, currentTag)
                    
                    currentTag = tag.lower()
                    relevant = []
                    print("tag changed to", currentTag)
                    
                relevant.append(line.split(" "))
        print(resultCode)

run()






