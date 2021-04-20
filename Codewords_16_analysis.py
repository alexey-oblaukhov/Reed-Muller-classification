
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
    
#Function which checks whether the set A has a 3-flat as its subset.
#When A is the codeword of weight 16, having a 3-flat as a subset is equivalent 
#to the statement "A can be represented as the union of 2 non-intersecting 3-flats".
#So in this case, any vector of A lies in some 3-flat, and we can start our search 
#from any vector of A
def does_have_3flat(A):
    
    output = []
    
    #We pick the first vector (arbitrary vector), and try to construct 
    #all possible 1-, 2- and then 3-flats which contain it.
    #If we succeed to find a 3-flat, we return 'True'
    v_b = A[0]
    flat = [v_b]
    
    for v_1 in A:
        if v_1 in flat:
            continue
            
        flat.append(v_1)
        diff_1 = v_1^v_b
        
        for v_2 in A:
            if v_2 in flat:
              continue
            if (v_2^diff_1) not in A:
              continue
            flat.append(v_2)
            flat.append(v_2^diff_1)
            diff_2 = v_2^v_b
            
            for v_3 in A:
                if v_3 in flat:
                    continue
                if ((v_3^diff_1) not in A) or ((v_3^diff_2) not in A) or ((v_3^diff_1^diff_2) not in A):
                    continue
                   
                #If we have reached this point, v_b + <v_1,v_2,v_3> form a 3-flat
                return True
            
            flat.pop()
            flat.pop()
        
        flat.pop()
    
    return False
        
#This function finds all possible representations
#of a codeword of weight 12 as a symm. difference of 2 3-flats
#For each such representation, there exist 6 parallel edges, 
#which can be arbitrarily split into two groups of 3 edges
#so that each group is contained in one of the two 3-flats
def decompose_12_into_2_3flats(A):
    
    if (len(A) != 12):
        print("Error, size is not 12!")
        return
    
    #First we find all possible edge directions
    directions = set()
    for u in A:
        for v in A:
            if (u != v):
                directions.add(u^v)
    
    decompositions = []
    #Then we search for groups of 6 parallel edges
    for dir in directions:
        pairs = []
        for v in A:
            if ((v^dir) in A) and (v < v^dir):
                pairs.append(v)
        
        if (len(pairs) == 6):
            #found 6 parallel edges, checking possible representations
            triples = list(itertools.combinations(pairs, 3))
            
            for triple in triples:
                v_a = 0
                v_b = 0
                for v in pairs:
                    if v in triple:
                        v_a ^= v
                    else:
                        v_b ^= v
                
                if (v_a == v_b) or (v_a == v_b^dir):
                    #Found a representation as a symmetric difference, adding it
                    decomposition = [[],[]]
                    for v in pairs:
                        if v in triple:
                            decomposition[0] += [v,v^dir]
                        else:
                            decomposition[1] += [v,v^dir]
                    
                    decomposition[0] += [v_a,v_a^dir]
                    decomposition[1] += [v_b,v_b^dir]
                    
                    decompositions.append([set(decomposition[0]),set(decomposition[1])])
                
                    
    return decompositions

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
#of a codeword A of weight 16 into a symmetric difference of 3 3-flats
#For any codeword of weight 16 which is a symmetric difference of 3 3-flats,
#for at least one of these 3-flats, exactly 6 vectors of it are in A
#so we can decompose A into a symmetric difference of a 3-flat and a codeword of weight 12
#by searching for groups of 3 parallel edges such that the fourth edge,
#which complements these 3 edges to a 3-flat, is not in A.
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
                #B is a codeword of weight 12 which we will decompose into 2 3-flats
                B = A.copy()
                B += [triple[0]^triple[1]^triple[2],triple[0]^triple[1]^triple[2]^dir]
                B.remove(triple[0]); B.remove(triple[1]); B.remove(triple[2])
                B.remove(triple[0]^dir); B.remove(triple[1]^dir); B.remove(triple[2]^dir)
                
                decomps_of_B = decompose_12_into_2_3flats(B)
                
                if (len(decomps_of_B) > 0):
                    #construct all possible decompositions into 3 3-flats, checking for repetitions
                    flat_1 = set([triple[0],triple[1],triple[2],triple[0]^triple[1]^triple[2],triple[0]^dir,triple[1]^dir,triple[2]^dir,triple[0]^triple[1]^triple[2]^dir])
                    for decomp in decomps_of_B:
                        new_decomp = [flat_1,decomp[0],decomp[1]]
                        is_good = True
                        for existing_decomp in decompositions:
                            if equal_ignore_order(existing_decomp,new_decomp):
                                is_good = False
                                break
                        
                        if is_good:
                            decompositions.append(new_decomp)
                
    return decompositions

#Function which counts the decompositions of the form 2 of the given codeword.
#It does so by checking sizes of pairwise intersections between A1, A2 and A3.
#Representations of the form 2 have 2 intersections of size 2 and 1 of size 0.
def count_representations_of_type_2(decompositions): 

    count = 0
    
    for triple in decompositions:
        intersection = [len(triple[0]&triple[1]),len(triple[0]&triple[2]),len(triple[1]&triple[2])]
        if ((intersection.count(2) == 2) and (intersection.count(0) == 1)):
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

#Construct representatives for each of the 12 equivalence classes 
c1 = symm_diff(A12, flat(e(4),e(3),e(1),e(2)))
c2 = symm_diff(A12, flat(e(4),e(3),e(1)^e(5),e(2)))
c3 = symm_diff(A12, flat(e(4),e(3),e(1),e(6)))
c4 = symm_diff(A12, flat(e(4),e(3),e(1)^e(5),e(6)))
c5 = symm_diff(A12, flat(e(4),e(3),e(6),e(7)))

c6 = symm_diff(A12, flat(e(4),e(5),e(1),e(2)))
c7 = symm_diff(A12, flat(e(4),e(5),e(1),e(6)))
c8 = symm_diff(A12, flat(e(4),e(5),e(6),e(7)))

c9 = symm_diff(A12, flat(e(4),e(1)^e(4),e(2)^e(5),e(6)))
c10 = symm_diff(A12, flat(e(4),e(1)^e(4),e(6),e(7)))

c11 = symm_diff(A12, flat(0,e(1),e(4),e(2)^e(5)))
c12 = symm_diff(A12, flat(0,e(1),e(4),e(6)))

#First we check which of the codewords contain a 3-flat
print("c1 " + (does_have_3flat(c1) and "contains" or "does not contain") + " a 3-flat")
print("c2 " + (does_have_3flat(c2) and "contains" or "does not contain") + " a 3-flat")
print("c3 " + (does_have_3flat(c3) and "contains" or "does not contain") + " a 3-flat")
print("c4 " + (does_have_3flat(c4) and "contains" or "does not contain") + " a 3-flat")
print("c5 " + (does_have_3flat(c5) and "contains" or "does not contain") + " a 3-flat")
print("c6 " + (does_have_3flat(c6) and "contains" or "does not contain") + " a 3-flat")
print("c7 " + (does_have_3flat(c7) and "contains" or "does not contain") + " a 3-flat")
print("c8 " + (does_have_3flat(c8) and "contains" or "does not contain") + " a 3-flat")
print("c9 " + (does_have_3flat(c9) and "contains" or "does not contain") + " a 3-flat")
print("c10 " + (does_have_3flat(c10) and "contains" or "does not contain") + " a 3-flat")
print("c11 " + (does_have_3flat(c11) and "contains" or "does not contain") + " a 3-flat")
print("c12 " + (does_have_3flat(c12) and "contains" or "does not contain") + " a 3-flat")
print("-----------------------------------------")

#c1-c4 contain a 3-flat, so we do not consider them further,
#as they have representation in the form 1 (two non-intersecting 3-flats).
#Next we find all distinct decompositions of remaining codewords
#into symmetric difference of 3 3-flats
c5_dec = decompose_into_3_3flats(c5)
c6_dec = decompose_into_3_3flats(c6)
c7_dec = decompose_into_3_3flats(c7)
c8_dec = decompose_into_3_3flats(c8)
c9_dec = decompose_into_3_3flats(c9)
c10_dec = decompose_into_3_3flats(c10)
c11_dec = decompose_into_3_3flats(c11)
c12_dec = decompose_into_3_3flats(c12)

#Now we check that the codewords c9-c12 have decompositions in the form 2
print("c9 " + ((count_representations_of_type_2(c9_dec) > 0) and "has" or "does not have") + " representations of the form 2")
print("c10 " + ((count_representations_of_type_2(c10_dec) > 0) and "has" or "does not have") + " representations of the form 2")
print("c11 " + ((count_representations_of_type_2(c11_dec) > 0) and "has" or "does not have") + " representations of the form 2")
print("c12 " + ((count_representations_of_type_2(c12_dec) > 0) and "has" or "does not have") + " representations of the form 2")
print("-----------------------------------------")

#For the remaining four codewords (of classes 5-8),
#we check that they all have distinct number of decompositions
#into a 3 3-flats, so they are inequivalent.
#We also count the number of decompositions of the form 2,
#since that will be used when enumerating and counting the codewords. 
print("c5 has " + str(len(c5_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_2(c5_dec)) + " of which are of the form 2")
print("c6 has " + str(len(c6_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_2(c6_dec)) + " of which are of the form 2")
print("c7 has " + str(len(c7_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_2(c7_dec)) + " of which are of the form 2")
print("c8 has " + str(len(c8_dec)) + " distinct unordered decompositions into 3 3-flats, " + str(count_representations_of_type_2(c8_dec)) + " of which are of the form 2")
print("-----------------------------------------")

print("Press 'Enter' to exit")
input()
