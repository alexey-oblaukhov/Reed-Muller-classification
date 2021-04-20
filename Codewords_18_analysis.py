
import itertools

#Function which generates a basis vector
def e(i):
    return 2**i
    
#Function which generates a 3-flat with shift 'a' and basis {b1,b2,b3}
def flat(a,b1,b2,b3):
    return [a ^ v for v in [0,b1,b2,b3,b1^b2,b1^b3,b2^b3,b1^b2^b3]]
   
#Function for the symmetric difference of lists A1, A2   
def symm_diff(A1,A2):
    output = A1.copy()
    for v in A2:
        if v in output:
            output.remove(v)
        else:
            output.append(v)
    
    return output

#This function finds the representation (it is unique)
#of a codeword of weight 14 as a symm. difference of 2 3-flats.
#The algorithm finds a group of 3 parallel edges, 
#(which must belong to one of the 3-flats)
#then finds a vector which is the sum of some 3 vectors from these edges' ends,
#(which is the 7th vector of the 3-flat)
#then sums the remaining three vectors to obtain the missing 8th vector
def decompose_14_into_2_3flats(A):

    if (len(A) != 14):
        print("Error, size is not 14!")
        return
    
    
    #First we find all possible edge directions
    directions = set()
    for u in A:
        for v in A:
            if (u != v):
                directions.add(u^v)
    
    #Then we search for groups of 3 parallel edges
    for dir in directions:
        pairs = []
        for v in A:
            if ((v^dir) in A) and (v < v^dir):
                pairs.append(v)
        
        if (len(pairs) == 3):
            #If found 3 parallel edges, decompose further
            
            decomposition = [[],[]]
            decomposition[0] += [pairs[0],pairs[1],pairs[2],pairs[0]^dir,pairs[1]^dir,pairs[2]^dir]

            for u in A:
                if ((pairs[0]^pairs[1]^pairs[2] == u) or (pairs[0]^pairs[1]^pairs[2]^dir == u)):
                    decomposition[0] += [u,u^dir] 
                    decomposition[1] = symm_diff(A, decomposition[0])

                    return [set(decomposition[0]),set(decomposition[1])]

#This function is used to check whether two decompositions
#into a symmetric difference of 3 3-flats are the same
def equal_ignore_order(a, b):
    unmatched = list(b)
    for element in a:
        try:
            unmatched.remove(element)
        except ValueError:
            return False
    return not unmatched 

#The function which searches for all distinct decompositions 
#of a codeword A of weight 18 into a symmetric difference of 3 3-flats
#For any codeword of weight 18 which is a symmetric difference of 3 3-flats,
#for at least one of these 3-flats, exactly 6 vectors of it are in A
#so we can decompose A into a symmetric difference of a 3-flat and a codeword of weight 14
#by searching for groups of 3 parallel edges such that the fourth edge,
#which complements these 3 edges to a 3-flat, is not in A
def decompose_into_3_3flats(A):

    decompositions = []
    directions = set()
    
    for u in A:
        for v in A:
            if (u != v):
                directions.add(u^v)
    
    #this dict will contain all groups of >=3 parallel edges
    groups_of_parallel_edges = {}
    for dir in directions:
        pairs = []
        for v in A:
            if ((v^dir) in A) and (v < v^dir):
                pairs.append(v)
        
        if (len(pairs) >= 3):
            groups_of_parallel_edges[dir] = pairs

    for dir in groups_of_parallel_edges:
        #check all groups of 3 parallel edges
        triples = list(itertools.combinations(groups_of_parallel_edges[dir], 3))
        
        for triple in triples:
            if not (triple[0]^triple[1]^triple[2] in A) and not (triple[0]^triple[1]^triple[2]^dir in A):
                #if the complement of these 3 edges to a 3-flat is not in A, we have a decomposition
                #B is a codeword of weight 14 which we will decompose into 2 3-flats
                B = A.copy()
                B += [triple[0]^triple[1]^triple[2],triple[0]^triple[1]^triple[2]^dir]
                B.remove(triple[0]); B.remove(triple[1]); B.remove(triple[2])
                B.remove(triple[0]^dir); B.remove(triple[1]^dir); B.remove(triple[2]^dir)
                
                decomp_of_B = decompose_14_into_2_3flats(B)
                
                flat_1 = set([triple[0],triple[1],triple[2],triple[0]^triple[1]^triple[2],triple[0]^dir,triple[1]^dir,triple[2]^dir,triple[0]^triple[1]^triple[2]^dir])
                
                new_decomp = [flat_1,decomp_of_B[0],decomp_of_B[1]]
                is_good = True
                for existing_decomp in decompositions:
                    if equal_ignore_order(existing_decomp,new_decomp):
                        is_good = False
                        break
                    
                if is_good:
                    decompositions.append(new_decomp)
    
    return decompositions
    

#Function which counts the decompositions of the form 1 of the given codeword.
#It does so by checking sizes of pairwise intersections between A1, A2 and A3.
#Representations of the form 1 have 2 intersections of size 2 and 1 of size 1.
def count_representations_of_type_1(decompositions): 

    count = 0
    
    for triple in decompositions:
        intersection = [len(triple[0]&triple[1]),len(triple[0]&triple[2]),len(triple[1]&triple[2])]
        if ((intersection.count(2) == 2) and (intersection.count(1) == 1)):
            count += 1
        
    return count


#########################
####### MAIN BODY #######
#########################

#Many obtained equivalence classes have the same A1 and A2 3-flats,
#So we put their symmetric difference into A12 
A1 = flat(0,e(1),e(2),e(3))
A2 = flat(0,e(3),e(4),e(5))
A12 = symm_diff(A1,A2)

#Construct representatives for each of the 10 equivalence classes 
c1 = symm_diff(A12, flat(0,e(4),e(5)^e(1),e(6)))
c2 = symm_diff(A12, flat(0,e(4),e(6),e(7)))

c3 = symm_diff(A12, flat(e(4),e(1),e(2),e(6)))
c4 = symm_diff(A12, flat(e(4),e(1)^e(5),e(2),e(6)))
c5 = symm_diff(A12, flat(e(4),e(1),e(6),e(7)))
c6 = symm_diff(A12, flat(e(4),e(1)^e(5),e(6),e(7)))
c7 = symm_diff(A12, flat(e(4),e(6),e(7),e(8)))

A2 = flat(0,e(4),e(5),e(6))
A12 = symm_diff(A1,A2)

c8 = symm_diff(A12, flat(e(1),e(1)^e(4),e(2)^e(5),e(3)^e(6)))
c9 = symm_diff(A12, flat(e(1),e(1)^e(4),e(2)^e(5),e(7)))
c10 = symm_diff(A12, flat(e(1),e(1)^e(4),e(7),e(8)))

#First we construct all distinct decompositions of the codewords
#into symmetric difference of 3 3-flats
c1_dec = decompose_into_3_3flats(c1)
c2_dec = decompose_into_3_3flats(c2)
c3_dec = decompose_into_3_3flats(c3)
c4_dec = decompose_into_3_3flats(c4)
c5_dec = decompose_into_3_3flats(c5)
c6_dec = decompose_into_3_3flats(c6)
c7_dec = decompose_into_3_3flats(c7)
c8_dec = decompose_into_3_3flats(c8)
c9_dec = decompose_into_3_3flats(c9)
c10_dec = decompose_into_3_3flats(c10)

#Check that the codewords c3-c6 have decompositions of the the form 1
print("c3 " + ((count_representations_of_type_1(c3_dec) > 0) and "has" or "does not have") + " representations of the form 1")
print("c4 " + ((count_representations_of_type_1(c4_dec) > 0) and "has" or "does not have") + " representations of the form 1")
print("c5 " + ((count_representations_of_type_1(c5_dec) > 0) and "has" or "does not have") + " representations of the form 1")
print("c6 " + ((count_representations_of_type_1(c6_dec) > 0) and "has" or "does not have") + " representations of the form 1")
print("-----------------------------------------")

#for the remaining codewords of classes 1,2,7,8,9,10
#check that they all have distinct number of decompositions
#into three 3-flats, so they are inequivalent
print("c1 has " + str(len(c1_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_1(c1_dec)) + " of which are of the form 1")
print("c2 has " + str(len(c2_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_1(c2_dec)) + " of which are of the form 1")
print("c7 has " + str(len(c7_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_1(c7_dec)) + " of which are of the form 1")
print("c8 has " + str(len(c8_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_1(c8_dec)) + " of which are of the form 1")
print("c9 has " + str(len(c9_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_1(c9_dec)) + " of which are of the form 1")
print("c10 has " + str(len(c10_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_1(c10_dec)) + " of which are of the form 1")
print("-----------------------------------------")

print("Press 'Enter' to exit")
input()
