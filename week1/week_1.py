import sys                                             #for argv[1]

def line_to_tokens(line) :                             # function for getting elemts , nodes, and value
    
    words = line.split()
    element_name = words[0]                            #element
    node_1 = words[1]                                  #node
    node_2 = words[2]                                  #node
    value = words[3]                                   #value
    
    ans = [element_name, node_1, node_2, value]        #temporary list
    return ans

def reversing(reverse):
    for x in reverse[::-1]:                            #reversing the list
        for y in x[::-1]:                              #reversing list elements
            print(y, end=' ')
        print('')
    print('')

#print(len(sys.argv))                                   #testing

try :
 len(sys.argv) == 2
 
 file = sys.argv[1]                                     #passing the file for easy usage
 
 if(not file.endswith(".netlist")) :                    #checking if netlist file has been given
         print("Invalid file type")
         
 else :
     containorlist =[]                                   #dummy list for storing lines one by one
     with open(file, "r") as F :                         #using with takes care of closing the file
         for lines in F.readlines():
            
             containorlist.append(lines.split('\n')[0])   #getting lines here splitting results in 2 set list with '' as other element so using [0] to remove it
         
         #print(containorlist)                            #checking
        
         try :                                            #using try to ensure there is .circuit in list
             circuit_position = containorlist.index(".circuit")
             end_position = containorlist.index(".end")

             required_list = containorlist[circuit_position+1:end_position]

             finalcontainer=[]                             #dummy container to store final answer from analysing tokenwise

             for line in required_list : 
             
                finalcontainer.append(line_to_tokens(line))

             #print(finalcontainer)                         #testing

             #print(finalcontainer)                         #testing

             reversing(finalcontainer)                      #reverse printing 
         
         except Exception :                                 #when second argument : netlist is absent
             print("Invalid file arguments, .circuit or .end is absent")
        
 
      

except  Exception : 
        print("Invalid command line, please enter a netlist file")

