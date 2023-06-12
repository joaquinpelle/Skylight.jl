 # Steps to translate from C:
 1. Change function to `metric!` and arguments to `g`, `position`, `spacetime::SuperposedPNSpacetime`
 2. Extract `t,x,y,z = position`
 3. Remove all variable declarations
 4. Remove all remaining "double"
 5. Remove all ";" (be careful there are no two or more commands in one line)
 6. Replace all // by  (checkout for  C macros before)
 7. Replace `result[i]` by `resulti`
 8. Use regular expressions to replace `pow(base, exp)` as follows (ChatGPT is very helpful to produce and interpret these):
    8.1. Find pow\(([^,]*),\s*([^)]+)\) using Reg Ex
    8.2. Replace by ($1)^($2) 
 9. Replace gcov by g
 10. Replace `return"` by `return nothing"`
 11. Remove all unused variables as indicated by the IDE
 12. Create setup function for variables imported from bbh_params struct
 13. In case of keeping arrays, shift indexes because Julia is 1-based (C is 0-based)
     To do this you can search for the regular expression (\w+)\[\s*(\d+)\s*\] and replace by $1[$2+1]
     Unfortunately regex's cannot do arithmetic operations, but the compiler will take care of them
     In case there is [][] or more general indexing you need to generalize the regex
 14. Remove final } and add end
 15. Replace "++" by "+=1"
 16. Replace integers that are used in ifs and whiles by booleans (or minimally you can change by e.g. while(keep_iterating==1))
