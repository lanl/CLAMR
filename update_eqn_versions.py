import itertools

""" This program spits out many of the parenthesizations of 

    U_new = U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus)
                     - wminusx_H + wplusx_H - wminusy_H + wplusy_H;

    U_new = U - (deltaT / dr) * (F_plus - F_minus + G_plus - G_minus) - wminusx_H + wplusx_H - wminusy_H + wplusx_H;


    Particularly, it makes all the parenthesizations of (F_plus - F_minus + G_plus - G_minus) and combines them with all
    possible parenthesizations of (- wminusx_H + wplusx_H  - wminusy_H + wplusy_H).  The parenthesization of 
    (F_plus - F_minus + G_plus - G_minus) is multiplied by (deltaT/dr) and then subtracted from U.  And then finally, the 
    parenthesization of the sum (-wminusx_H + wplus_h - wminusy_H + wplus_H) is added to this result. 

    The function allbinarytrees was taken from 

    https://stackoverflow.com/questions/17242117/what-python-code-generates-all-possible-groupings-trees-for-binary-operators

    Written by Vanessa Job, August 2019.

    UPDATE: Now the program generates all possible orderings for F_plus - F_minus + G_plus - G_minus and then does all possible 
    parenthesizations of this AND then generates all possible orderings of - wminusx_H + wplusx_H - wminusy_H + wplusx_H with 
    all possible parenthesizations.   These are combined to make equations and the resulting equations are output. """ 

def all_permuted_equations(variables): 
    """ Given a list of single character variable names, return permuations of the sums of these variables.  These are 
        constricte to single character varaible names because if not, it will break allbinarytrees.  V, to fix allbinarytrees.   """
    all_perms = list(itertools.permutations(variables))
    all_equations = list()
    for thing in all_perms:
        new_eq = ""
        for i in range(len(thing) - 1): 
            new_eq += thing[i] + "+"
        new_eq += thing[-1] 
        all_equations.append(new_eq)

    return all_equations    
             
    
def allbinarytrees(s):
    """ Return all possible parenthesizations of the expression s which is assumed to be totally without parentheses. """
    if len(s) == 1:
        yield s
    else:
        for i in range(1, len(s), 2):
            for l in allbinarytrees(s[:i]):
                for r in allbinarytrees(s[i+1:]):
                    yield '({}{}{})'.format(l, s[i], r)



if __name__ == "__main__":
 
    all_perms = all_permuted_equations(["a", "b", "c", "d"])
    #print(all_perms)

    #To see that we have all the possible parenthesizations:
    #
    #for equation in all_perms:
    #     for t in allbinarytrees(equation):
    #        print(t)
    #     print("\n")


    #print("All equations: ")

    eqn_file = open("update_eqn_versions.h","w")

    eqn_file.write("inline real_t U_fullstep_version(\n")
    eqn_file.write("        real_t    deltaT,\n")
    eqn_file.write("        real_t    dr,\n")
    eqn_file.write("        real_t    U,\n")
    eqn_file.write("        real_t    F_plus,\n")
    eqn_file.write("        real_t    F_minus,\n")
    eqn_file.write("        real_t    G_plus,\n")
    eqn_file.write("        real_t    G_minus,\n")
    eqn_file.write("        real_t    wplusx_H,\n")
    eqn_file.write("        real_t    wminusx_H,\n")
    eqn_file.write("        real_t    wplusy_H,\n")
    eqn_file.write("        real_t    wminusy_H\n")
    eqn_file.write(") {")

    count = 0
    all_expression = list() 
    for perm_1 in all_perms:  
        for perm_2 in all_perms:
            for first_part in allbinarytrees(perm_1):
                expression_1 = str(first_part)
                expression_1 = expression_1.replace('a', ' F_plus ' )
                expression_1 = expression_1.replace('b', ' -F_minus '  )
                expression_1 = expression_1.replace('c','  G_plus '  )
                expression_1 = expression_1.replace('d', '-G_minus '  )
        
                for second_part in allbinarytrees(perm_2):
                    expression_2 = str(second_part)
                    expression_2 = expression_2.replace('a', ' -wminusx_H '  ) 
                    expression_2 = expression_2.replace('b', ' wplusx_H ')
                    expression_2 = expression_2.replace('c', ' -wminusy_H ')
                    expression_2 = expression_2.replace('d', ' wplusy_H ')
        
                    eqn_file.write("\n")
                    ifstr = "#if ( UPDATE_EQUATION_VERSION ==  " + str(count) + " )\n"
                    eqn_file.write(ifstr)
                    expression = "return  U -  (deltaT/dr) * (" + expression_1 + ") + (" +  expression_2 + ");\n"
                    eqn_file.write(expression)
                    eqn_file.write("#endif\n")
                    count += 1
    
    eqn_file.write("//If UPDATE_EQUATION_VERSION was not equal to any listed options, return default equation version\n")
    eqn_file.write("return U - (deltaT / dr)*(F_plus - F_minus + G_plus - G_minus) + - wminusx_H + wplusx_H  + - wminusy_H + wplusy_H;\n")
    eqn_file.write("}")
