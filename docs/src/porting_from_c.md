 # Hints to translate a metric from C:
 1. Change function to `metric!` and arguments to `g`, `position`, `spacetime::SPNOldSpacetime`
 2. Extract `t,x,y,z = position`
 3. Remove all variable declarations
 4. Remove all remaining "double"
 5. Remove all ";" (be careful there are no two or more commands in one line)
 6. Replace all // by  (checkout for  C macros before)
 7. Replace `result[i]` by `resulti`
 8. Use regular expressions to replace `pow(base, exp)` as follows (ChatGPT is very helpful to produce and interpret these):
    8.1. Find pow\(([^,]*),\s*([^)]+)\) using Reg Ex
    8.2. Replace by ($1)^($2) 
 9. Replace gcov by g (just for consistency)
 10. Replace `return"` by `return nothing"`
 11. Remove all unused variables as shown by the IDE
 12. Create setup function for variables imported from bbh_params struct
 13. In case of keeping arrays, shift indexes because Julia is 1-based (whereas C is 0-based)
     To do this you can search for the regular expression (\w+)\[\s*(\d+)\s*\] and replace by $1[$2+1]
     After that you can search all `(1+1)` and replace by `2` and so on.
     In case there is `[][]` or other type of indexing you need to generalize the regex
 14. Remove final `}` and add `end`
 15. Replace `++` by `+=1`
 16. Replace integers that are used in `if`s and `while`s by booleans (or minimally you can change by e.g. `while(keep_iterating==1))` and then replace `keep_iterating=1` by `keep_iterating=true` and `keep_iterating=0` by `keep_iterating=false`
