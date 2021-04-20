
#i_XYZ stands for the size of the intersection of the 3-flats A_X, A_Y and A_Z
#Any pairwise intersection can only have cardinality at most 2
#Any triple intersection cannot be greater than all its pairwise sub-intersections
#WLOG (renumbering the sets if we need to) we can assume that
#i_12 is the largest pairwise intersection, and also
#i_12 >= i_13 >= i_14

print("'Intersection configurations' of four 3-flats:")

counter = 0

for i_12 in range(3):
    for i_13 in range(i_12 + 1):
        for i_14 in range(i_13 + 1):
            for i_23 in range(i_12 + 1):
                for i_24 in range(i_12 + 1):
                    for i_34 in range(i_12 + 1):
                        for i_123 in range(min(i_12,i_13,i_23) + 1):
                            for i_124 in range(min(i_12,i_14,i_24) + 1):
                                for i_134 in range(min(i_13,i_14,i_34) + 1):
                                    for i_234 in range(min(i_23,i_24,i_34) + 1):
                                        for i_1234 in range(min(i_123,i_124,i_134,i_234) + 1):
                                            
                                            #Calculating the weight of the codeword A constructed as the symmetric difference of four 3-flats intersecting in such ways
                                            wt = 32 - 2 * (i_12 + i_13 + i_14 + i_23 + i_24 + i_34) + 4 * (i_123 + i_124 + i_134 + i_234) - 8 * i_1234
                                            
                                            if (wt == 18):
                                                #Calculating the weight of the intersection of each of the 3-flats with the resulting codeword A
                                                w_1 = 8 - (i_12 + i_13 + i_14) + 2 * (i_123 + i_124 + i_134) - 4 * i_1234
                                                w_2 = 8 - (i_12 + i_23 + i_24) + 2 * (i_123 + i_124 + i_234) - 4 * i_1234
                                                w_3 = 8 - (i_13 + i_23 + i_34) + 2 * (i_123 + i_134 + i_234) - 4 * i_1234
                                                w_4 = 8 - (i_14 + i_24 + i_34) + 2 * (i_124 + i_134 + i_234) - 4 * i_1234
                                                
                                                #If some 3-flat already has >= 6 vectors in A, we do not need to consider such case
                                                if (w_1 < 6) and (w_2 < 6) and (w_3 < 6) and (w_4 < 6):
                                                    
                                                    #Next we check if for some 3-flat other than A_1, its 3 pairwise intersections with other sets
                                                    #"majorize" the intersections of A_1 with the other 3 sets. If that's the case, we don't consider
                                                    #this case since it is a repetition of another already considered case
                                                    
                                                    alt = [i_12,i_23,i_24]
                                                    alt.sort(reverse = True)
                                                    if (alt[0] >= i_12) and (alt[1] >= i_13) and (alt[2] >= i_14) and ((alt[0] > i_12) or (alt[1] > i_13) or (alt[2] > i_14)):
                                                        #majorized by A_2, ignore
                                                        continue
                                                    
                                                    alt = [i_13,i_23,i_34]
                                                    alt.sort(reverse = True)
                                                    if (alt[0] >= i_12) and (alt[1] >= i_13) and (alt[2] >= i_14) and ((alt[0] > i_12) or (alt[1] > i_13) or (alt[2] > i_14)):
                                                        #majorized by A_3, ignore
                                                        continue
                                                        
                                                    alt = [i_14,i_24,i_34]
                                                    alt.sort(reverse = True)
                                                    if (alt[0] >= i_12) and (alt[1] >= i_13) and (alt[2] >= i_14) and ((alt[0] > i_12) or (alt[1] > i_13) or (alt[2] > i_14)):
                                                        #majorized by A_4, ignore
                                                        continue
                                                    
                                                    counter += 1
                                                    print(str(counter) + ". " + str([i_12,i_13,i_14,i_23,i_24,i_34,i_123,i_124,i_134,i_234,i_1234]))

#The program outputs 10 configurations.
#However, it is easy to see that some of the configurations are the same
#if we renumber the sets. More precisely:
#Configurations 1 and 2 are the same (switch A_3 and A_4)
#Configurations 4 and 5 are the same (switch A_2 and A_3)
#Configurations 8 to 10 are the same (switch A_2, A_3 and A_4 around)
#Configurations 1, 7 and 8 are the same, (switch flats 1-4 and then 2-3)
#It is easy to see that all remaining configurations are not the same.
#So unique ones are (for example): 2, 3, 5, 6, 7, 10.

print("Press 'Enter' to exit")
input()