import math
import cmath
import sys
import numpy as np

PI=np.pi  
CIRCUIT = ".circuit" 
END = ".end"
                                  
values_letters = {"T":1e+12,"t":1e+12,"G":1e+09 ,"g":1e+09 ,"Meg":1e+06 ,"meg": 1e+06,           #Dictionary to get multiplier 
                  "K":1e+03 ,"k":1e+03, "M":1e-03 ,"m":1e-03, "U":1e-06 ,"µ": 1e-06,             #values are eneterd in kilo, micro
                  "u":1e-06, "n":1e-09,"N": 1e-09, "p":1e-12, "P": 1e-12, "f":1e-15 , "F":1e-15} #,etc


def isnum(var):               
#This function is used to verify wheter value for each element is in correct format or not
    try :
        var = float(var)
    except Exception :return False
    return not math.isnan (var)

def checkRepetition(list):
#Checks if same element names occur twice. NOTE: The program will issue a warning
#but continue operation assuming the elements are different
    checking_list = []
    for obj in list:
        if obj.name not in checking_list : checking_list.append(obj.name)
        else :
            print("----------\n")
            print("Warning ! : Element "+obj.name+" is defined twice in netlist. We will assume them to be distinct and solve. If not intended, correct this mistake")
            print("----------\n")

# Defining math function for ease of ease
def sin(val):                            
    return math.sin(val)
def cos(val):
    return math.cos(val)

#Defining Eulers form of complex numbers
def polar(z):                           
    a=z.real
    b=z.imag
    if(math.hypot(a,b)>0.0001 or math.hypot(a,b)<(1.0e-10)):
      r = round(math.hypot(a,b),4)
      a=round(a,4)
      b=round(b,4)
    else :
      r= math.hypot(a,b)
    theta = "{:.2e}".format((math.atan2(b,a)) *180/PI)
        
    return r,theta, a,b

#Class component to store each element
class Components:                   
#Class componenet initialises each component depending on the element type. 
#For instance a DC voltage source wont have a phase value but AC would have phase intialised by classes constructor  
     def __init__(self,component_list):
        
        if len(component_list ) ==4:
            self.name = component_list[0]
            self.node1 = component_list[1]
            self.node2 = component_list[2]
            self.type = component_list[0][0]
            if(component_list[3][-1].isalpha()): 
                mult = values_letters[component_list[3][-1]]
                self.value = float(component_list[3][:-1]) * mult
            elif(isnum(component_list[3])) :self.value = float(component_list[3])
            else: print("ERROR : Invalid value encountered for element "+self.name);exit()
         
        elif len(component_list)==6 and component_list[3] == "ac":
            self.name = component_list[0]
            self.node1 = component_list[1]
            self.node2 = component_list[2]
            self.type = component_list[0][0]
            if(component_list[4][-1].isalpha()): 
                mult = values_letters[component_list[4][-1]]
                self.value = float(component_list[4][:-1]) * mult/2
            elif(isnum(component_list[4])) :self.value = float(component_list[4])/2
            else : print("ERROR : Invalid value encountered for element "+self.name);exit()
            self.phase = float(component_list[5])

           

        elif len(component_list) ==5 and component_list[3] == "dc":
            self.name = component_list[0]
            self.node1= component_list[1]
            self.node2 = component_list[2]
            if(component_list[4][-1].isalpha()): 
                mult = values_letters[component_list[4][-1]]
                self.value = float(component_list[4][:-1]) * mult
            elif(isnum(component_list[4])) :self.value = float(component_list[4])
            else : print("ERROR : Invalid value encountered for element "+self.name); exit()
            
            self.type = component_list[0][0]

        elif len(component_list) ==6 and (component_list[0][0] == "G" or component_list[0][0] == "E"):
            self.name = component_list[0]
            self.node1= component_list[1]
            self.node2 = component_list[2]
            self.VC_node1 = component_list[3]
            self.VC_node2 = component_list[4]
            if(component_list[5][-1].isalpha()): 
                mult = values_letters[component_list[5][-1]]
                self.value = float(component_list[5][:-1]) * mult
            elif(isnum(component_list[5])) :self.value = float(component_list[5])
            else : print("ERROR : Invalid value encountered for element "+self.name); exit()
            
            self.type = component_list[0][0]

        elif len(component_list) ==5 and (component_list[0][0] == "H" or component_list[0][0] == "F"):
            self.name = component_list[0]
            self.node1= component_list[1]
            self.node2 = component_list[2]
            self.CC_elem = component_list[3]
            if(component_list[4][-1].isalpha()): 
                mult = values_letters[component_list[4][-1]]
                self.value = float(component_list[4][:-1]) * mult
            elif(isnum(component_list[4])) :self.value = float(component_list[4])
            else : print("ERROR : Invalid value encountered for element "+self.name);exit()
            
            self.type = component_list[0][0]

     def print(self):
         print(self.name, end=" ")
         print(self.node1, end=" ")
         print(self.node2, end=" ")
         print(self.type, end=" ")
         print(self.value, end = " ")
#Get variable X vector
def get_x_vector(given_component_list):
#Here we are initailizing x vector of the eqn Mx =b               
    x_vector = []

    for elem in given_component_list:
        
        if ("V_"+elem.node1) not in x_vector:
            x_vector.append("V_"+elem.node1)

        if ("V_"+elem.node2) not in x_vector:
            x_vector.append("V_"+elem.node2)

        if elem.type == "V" or elem.type == "E" or elem.type =="H":
            x_vector.append("I_"+elem.name)

    return x_vector

#Get M and b for DC
#Algorithims used for the functions get_Matrix_AC and get_Matrix_dc are derived form lec10 of ee2015
def get_Matrix_dc(component_list, x_vector):                   
    len_of_matrix = len(x_vector)

    M_matrix = np.zeros((len_of_matrix, len_of_matrix), dtype =float)
    b_vector = np.zeros((len_of_matrix), dtype=float)

    for element in component_list:
 
        n1 = x_vector.index("V_"+element.node1)         
        n2 = x_vector.index("V_"+element.node2)

        if element.type == "R":
            M_matrix[n1, n1] += 1/element.value
            M_matrix[n2, n2] += 1/element.value
            M_matrix[n1, n2] -= 1/element.value
            M_matrix[n2, n1] -= 1/element.value

        elif element.type == "V":
            node_iv = x_vector.index("I_"+element.name)
            M_matrix[n1, node_iv] -= 1
            M_matrix[n2, node_iv] += 1
            M_matrix[node_iv, n1] += 1
            M_matrix[node_iv, n2] -= 1

            b_vector[node_iv] += element.value

        elif element.type == "I":
            b_vector[n1] -= element.value
            b_vector[n2] += element.value
        
        elif element.type == "L" :
            M_matrix[n1, n1] += 1e+100/element.value        #adding an extremely large value
            M_matrix[n2, n2] += 1e+100/element.value
            M_matrix[n1, n2] -= 1e+100/element.value
            M_matrix[n2, n1] -= 1e+100/element.value

        elif element.type == "C" :
            M_matrix[n1, n1] += 1e-100*element.value        #adding an extremely small value
            M_matrix[n2, n2] += 1e-100*element.value
            M_matrix[n1, n2] -= 1e-100*element.value
            M_matrix[n2, n1] -= 1e-100*element.value

        elif element.type == "G":
            n_m = x_vector.index("V_"+element.VC_node1) 
            n_n = x_vector.index("V_"+element.VC_node2)   
          
            M_matrix[n1,n_m] +=element.value
            M_matrix[n1,n_n] -=element.value
            M_matrix[n2,n_m] -=element.value
            M_matrix[n2,n_n] +=element.value

        elif element.type == "E":
            n_m = x_vector.index("V_"+element.VC_node1)
            n_n = x_vector.index("V_"+element.VC_node2)
            n_i = x_vector.index("I"+element.name)
           
            M_matrix[n1, n_i] +=1
            M_matrix[n2, n_i] -=1
            M_matrix[n_i, n1] +=1
            M_matrix[n_i, n2] -=1
            M_matrix[n_i, n_m] -= element.value
            M_matrix[n_i, n_n] += element.value
        
        elif element.type == "F":
            n_i = x_vector.index("I_V"+element.name)
            M_matrix[n1, n_i] += element.value
            M_matrix[n2, n_i] -= element.value

        elif element.type =="H":
            n_imn = x_vector.index("I_V"+element.name)
            n_ikl = x_vector.index("I_" + element.name)

            M_matrix[n1, n_ikl] +=1
            M_matrix[n2, n_ikl] -=1
            M_matrix[n_ikl, n1] +=1
            M_matrix[n_ikl, n2] -=1
            M_matrix[n_ikl, n_imn] -= element.value


    try:                                                     #used to check if GND node has been specified
         c=x_vector.index("V_GND")
         M_matrix[c,c]+=1
         return M_matrix, b_vector

    except Exception : 
        print("Error No GND FOUND! Please change 0 to GND, or recheck file")
        exit()
    
  

#Get M and b for AC
#Logic is similar to get_Matrix_ac except that for L and C we will take w into conductance calculation
def get_Matrix_AC(component_list, x_vector, frequency):       
    len_of_matrix = len(x_vector)

    M_matrix = np.zeros((len_of_matrix, len_of_matrix), dtype =complex)
    b_vector = np.zeros((len_of_matrix), dtype=complex)

    for element in component_list:

        n1 = x_vector.index("V_"+element.node1)
        n2 = x_vector.index("V_"+element.node2)

        if element.type == "R":
            M_matrix[n1, n1] += 1/element.value
            M_matrix[n2, n2] += 1/element.value
            M_matrix[n1, n2] -= 1/element.value
            M_matrix[n2, n1] -= 1/element.value

        if element.type == "L":
            M_matrix[n1, n1] += 1/(element.value*2*PI*frequency*(1j))
            M_matrix[n2, n2] += 1/(element.value*2*PI*frequency*(1j))
            M_matrix[n1, n2] -= 1/(element.value*2*PI*frequency*(1j))
            M_matrix[n2, n1] -= 1/(element.value*2*PI*frequency*(1j))

        if element.type =="C":
            M_matrix[n1, n1] += (element.value*2*PI*frequency*(1j))
            M_matrix[n2, n2] += (element.value*2*PI*frequency*(1j))
            M_matrix[n1, n2] -= (element.value*2*PI*frequency*(1j))
            M_matrix[n2, n1] -= (element.value*2*PI*frequency*(1j))
            

        if element.type == "V":
            node_iv = x_vector.index("I_"+element.name)
            M_matrix[n1, node_iv] -= 1
            M_matrix[n2, node_iv] += 1
            M_matrix[node_iv, n1] += 1
            M_matrix[node_iv, n2] -= 1
            try :
               b_vector[node_iv] +=( element.value * (cos(element.phase) + ( sin(element.phase) *(1j) ) ) )
            
            except Exception: b_vector[node_iv] = element.value

        if element.type == "I":
            try:
                b_vector[n1] -= (element.value * (cos(element.phase *PI/180) + ( sin(element.phase *PI/180) *(1j) ) ))
                b_vector[n2] += (element.value * (cos(element.phase *PI/180) + ( sin(element.phase *PI/180) *(1j) ) ))
            except Exception: 
                #element.print()
                b_vector[n1] -= (element.value )
                b_vector[n2] += (element.value ) 
            
         
        elif element.type == "G":
            n_m = x_vector.index("V_"+element.VC_node1) #m
            n_n = x_vector.index("V_"+element.VC_node2) #n
            #k=n1 l=n2
            M_matrix[n1,n_m] +=element.value
            M_matrix[n1,n_n] -=element.value
            M_matrix[n2,n_m] -=element.value
            M_matrix[n2,n_n] +=element.value

        elif element.type == "E":
            n_m = x_vector.index("V_"+element.VC_node1) #m
            n_n = x_vector.index("V_"+element.VC_node2) #n
            n_i = x_vector.index("I_"+element.name)#ivkl
            #k=n1 l=n2
            M_matrix[n1, n_i] +=1
            M_matrix[n2, n_i] -=1
            M_matrix[n_i, n1] +=1
            M_matrix[n_i, n2] -=1
            M_matrix[n_i, n_m] -= element.value
            M_matrix[n_i, n_n] += element.value
        
        elif element.type == "F":
            #n1=k, n2=l
            n_i = x_vector.index("I_V"+element.name)
            M_matrix[n1, n_i] += element.value
            M_matrix[n2, n_i] -= element.value

        elif element.type =="H":
            #n1 = k, n2=l
            n_imn = x_vector.index("I_V"+element.name)
            n_ikl = x_vector.index("I_" + element.name)

            M_matrix[n1, n_ikl] +=1
            M_matrix[n2, n_ikl] -=1
            M_matrix[n_ikl, n1] +=1
            M_matrix[n_ikl, n2] -=1
            M_matrix[n_ikl, n_imn] -= element.value
     
    try:
         c=x_vector.index("V_GND")
         M_matrix[c,c]+=1
         return M_matrix, b_vector

    except Exception : 
        print("Error No GND FOUND! Please change 0 to GND, or recheck file")
        exit()
    




# This function is used break down each line in of netlist into list containing that elements spcifics
def line_to_tokens(line):
    words = line.split()
    
    if len(words) == 4:                                                    #R,L,C
      element_name = words[0]                                              # element
      node_1 = words[1]                                                    # node
      node_2 = words[2]                                                    # node
      value = words[3]                                                     # value
      ans = [element_name, node_1, node_2, value]                          # temporary list
      return ans
    
    elif len(words) == 6 and ( (words[0][0] == "V") or words[0][0] == "I") :
       element_name = words[0]                                             # element
       node_1 = words[1]                                                   # node
       node_2 = words[2]                                                   # node
       type = words[3]                                                     #ac
       value = words[4]                                                    # value
       phase = words[5]                                                    #phase
       ans = [element_name, node_1, node_2, type ,value, phase]            # temporary list
       return ans

    elif len(words) == 5 and  ( (words[0][0] == "V") or words[0][0] == "I") :
    
       element_name = words[0]                                             # element
       node_1 = words[1]                                                   # node
       node_2 = words[2]                                                   # node
       type = words[3]                                                     #dc
       value = words[4]                                                    # value
       ans = [element_name, node_1, node_2, type ,value]                   # temporary list
       return ans      

    elif len(words) == 6 and  ( (words[0][0] == "E") or words[0][0] == "G"):
       element_name = words[0]                                             # element
       node_1 = words[1]                                                   # node
       node_2 = words[2]                                                   # node
       VS_node_1 = words[3]                                                #node corresponding controlling voltage
       VS_node_2 = words[4]                                                #node corresponding controlling voltage
       value = words[5]                                                    # value
       ans = [element_name,node_1, node_2, VS_node_1, VS_node_2 , value]   # temporary list
       return ans

    elif len(words) == 5 and  ( (words[0][0] == "H") or words[0][0] == "F"):
       element_name = words[0]                                             # element
       node_1 = words[1]                                                   # node
       node_2 = words[2]                                                   # node
       CC_element = words[3]                                               #controlling current
       value = words[4]                                                    # value
       ans = [element_name, node_1, node_2, CC_element ,value]             # temporary list
       return ans

    
#try to see if netlist has been passed    
try:
 len(sys.argv) == 2
 file = sys.argv[1]                                                        # passing the file for easy usage
except  Exception : 
        print("Invalid command line, please enter a netlist file")
        exit()

if(not file.endswith(".netlist")):                                         # checking if netlist file has been given
         print("Invalid file type, Please enter a .netlist file")
         exit()

else:
     containorlist = []                                                    # dummy list for storing lines one by one
     try:
         with open(file, "r") as F:                                        # using with takes care of closing the file
             for lines in F.readlines():
                 # getting lines here splitting results in 2 set list with '' as other element so using [0] to remove it
                 containorlist.append(lines.split("#")[0].split('\n')[0])
                 
             try:                                                          # using try to ensure there is .circuit and
                 circuit_position = containorlist.index(CIRCUIT)        # .end in list
                 end_position = containorlist.index(END)
             except Exception:                                             # when .circuit or .end is absent  is absent
                 print("Invalid file defination, .circuit or .end is absent")
         
             required_list = containorlist[circuit_position+1:end_position] #limiting the lines to contain those between .circuit and .en
             required_list= [i for i in required_list if i != '']           #This is used to remove any '' or blanks
             
             AC = False                                                     #Checking If AC or Not
             frequency =0
             for element in (containorlist):
                 if(element[:3] == ".ac") :
                     AC= True
                     
                     frequency = float(element.split()[-1])                 #storing frequency
    
             finalcontainer = []                                            # dummy container to store final answer from analysing
                                                                           # tokenwise
             for line in required_list:
                 finalcontainer.append(line_to_tokens(line))
            
             component_list = []                                            #List of objects of class components, it will have each
                                                                            #element in circuit stored in the form of object of 
                                                                            #class component for ease of use
             #Defining each object for each element present in circuit                                                              
             for a_component in finalcontainer:                             
                 def_object = Components(a_component)
                 component_list.append(def_object)
    
                 if def_object.type == "H" or def_object.type == "F":       #Here we try to check for CCXS
                     for elem in component_list:                            #As we will have to add a dummy voltage source of value 0
                         if elem.name == def_object.CC_elem:                #and introduce pseudo node 
    
                             index = component_list.index(elem)
                             component_list.remove(elem)
    
                             n1_elem = elem.node1
    
                             n_fake = "N"+ def_object.name
                             elem.node1 = n_fake
                             component_list.insert(index, elem)
                             #adding the dummy voltage
                             fakeVoltage = Components( [ ("V"+def_object.name), n1_elem, n_fake, "0"] )
                             component_list.append(fakeVoltage)
                            
             checkRepetition(component_list)                                 # This is used to check if any element name has repeated
    
             x_vector = get_x_vector(component_list)                         #Get variable list/vector
    
             if not AC:                                                      # For DC
                 M, b = get_Matrix_dc(component_list, x_vector)

                 try: 
                     x= np.linalg.solve(M,b)                                 #solves Mx =b and gives x
                 except np.linalg.LinAlgError :
                     print("The given circuit is insolvable. Please check the validity of given netlist")
                     exit()                                 
    
                 for i in range(len(x_vector)):                              #printing
    
                     if(x_vector[i][0] == "V"):
                         print("Voltage at node "+x_vector[i][( (x_vector[i]).index("_") +1):] + " : " , end=" ")
                         print(round(float(x[i]), 4), end=" V\n")
    
                    
                     elif(x_vector[i][0] == "I"):
                         print("Current through Voltage source " + x_vector[i][( (x_vector[i]).index("_") +1):] + " : ", end=" ")
                         print(round(float(x[i]), 4), end=" A\n")
    
             else :                                                         #For AC    
                 
                 M,b = get_Matrix_AC(component_list, x_vector,frequency)

                 try: 
                     x= np.linalg.solve(M,b)                                 #solves Mx =b and gives x
                 except np.linalg.LinAlgError :
                     print("The given circuit is insolvable. Please check the validity of given netlist")
                     exit()  
                 
                 
                 for i in range(len(x_vector)):                             #printing in cosine form
                     magnitude, phase, realPart, imaginaryPart = polar(x[i])
                     if(x_vector[i][0] == "V"):
                        
                         
                         print("Voltage at node "+x_vector[i][( (x_vector[i]).index("_") +1):] + " : " , end=" ")
    
                         if(magnitude==0): 
                             print( 0.0)
                         else :
                             print(str(magnitude) + "cos( " + "{:.2e}".format((2*PI*frequency))+"t+(" +(str(phase)) +")deg) V")
                             '''print("In complex form- ", end=" ")
                             print(complex(realPart,imaginaryPart))
                             print()'''
    
                     elif(x_vector[i][0] == "I"):
                         print("Current through Voltage source " + x_vector[i][( (x_vector[i]).index("_") +1):] + " : ", end=" ")
                         
                         if(magnitude==0) :print(0.0)
                         else:
                             print(str(magnitude) + "cos( " + "{:.2e}".format((2*PI*frequency))+"t+(" +(str(phase)) +")deg) A") 
                             '''print("In complex form- ", end=" ")
                             print(complex(realPart,imaginaryPart))      
                             print()    '''                                                                                    
                 print("\nIn Polar form\n")
    
                 for i in range(len(x_vector)):                                #printing in complex form
                     magnitude, phase, realPart, imaginaryPart = polar(x[i])
                     if(x_vector[i][0] == "V"):
                         print("Voltage at node "+x_vector[i][( (x_vector[i]).index("_") +1):] + " : " , end=" ")
    
                         if(magnitude==0): 
                             print( 0.0 )
                         else :
                             print(complex(realPart,imaginaryPart))
    
                     elif(x_vector[i][0] == "I"):
                         print("Current through Voltage source " + x_vector[i][( (x_vector[i]).index("_") +1):] + " : ", end=" ")
                         
                         if(magnitude==0) :print(0.0)
                         else:
                             print(complex(realPart,imaginaryPart))      
                                                                      
         
     except IOError:
         print("File not Found! Please recheck location")  