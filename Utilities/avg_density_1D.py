#A short script to create extended domains from the computational domains when there are 
#reflective boundary conditions.  Extended domains are not necessary for the calculations
#or publication, but sometimes are helpful for visualizing the system under investigation.
#Written by Daniel Frink, April 2020.


#reflect y-pos y-neg x-pos for dft_dens.dat
#reflect x_pos z_pos z_neg for dft_dens1.dat

file_name = input("Input file name: ")
out_file_name = input("Output file name: ")
num_dims = int(input("Enter number of dimensions: "))
startx = double(input("x value to start average: "))
avg_range = double(input("range of average: "))
Num_densities = double(input("number of density variables in file" ))

increment = 0.0
x_pos_bound = 0.0
y_pos_bound = 0.0
z_pos_bound = 0.0
x_neg_bound = 0.0
y_neg_bound = 0.0
z_neg_bound = 0.0

def storeData(file_name):
    output = []
    input_data = open(file_name,"r")
    for line in input_data:
        line = line.split()
        if(len(line)==1):
            continue
        else:
            output.append(line)
    #trim whitespace of end of data
    lastInd = len(output)-1
    while(len(output[lastInd])==0):
        output.pop(lastInd)
        lastInd-=1
    input_data.close()
    return output

def findAverage(startx,avg_range,Num_densities,data):

def reflect(axis,direction,data):
    global y_pos_bound
    global x_pos_bound
    global z_pos_bound
    global x_neg_bound
    global y_neg_bound
    global z_neg_bound
    global increment
    new_data = []
    if(axis == 0):
        num_chunks = (((y_pos_bound-y_neg_bound)/increment)+1)*(((z_pos_bound-z_neg_bound)/increment)+1)
        chunk_size = ((x_pos_bound-x_neg_bound)/increment)
        if(direction == 0):
            for line in data:
                if(line==[]):
                    continue
                reflected_line = line.copy()
                reflected_line[0] = float(reflected_line[0])
                if(reflected_line[0]==x_neg_bound):
                    continue
                else:
                    reflected_line[0] = -reflected_line[0]
                    reflected_line[0] = format(reflected_line[0],'.6f')
                    new_data.append(reflected_line)
            for i in range(0,int(num_chunks)):
                for j in range(0,int(chunk_size)):
                    data.insert((2*i*(int(chunk_size)+1)),new_data[(i*int(chunk_size))+j])
        elif(direction == 1):
            for line in data:
                if(line==[]):
                    continue
                reflected_line = line.copy()
                reflected_line[0] = float(reflected_line[0])
                if(reflected_line[0]==x_pos_bound):
                    continue
                else:
                    reflected_line[0] = x_pos_bound + (x_pos_bound-reflected_line[0])
                    reflected_line[0] = format(reflected_line[0],'.6f')
                    new_data.append(reflected_line)
            for i in range(0,int(num_chunks)):
                for j in range(0,int(chunk_size)):
                    data.insert((2*i*(int(chunk_size)+1)+int(chunk_size)+1),new_data[(i*int(chunk_size))+j])
    elif(axis==1):
        num_chunks = int((((z_pos_bound-z_neg_bound)/increment)+1)*((y_pos_bound-y_neg_bound)/increment)) #number of individual blocks
        y_chunks = int(((y_pos_bound-y_neg_bound)/increment))
        chunk_size = int(((x_pos_bound-x_neg_bound)/increment)+2) #size of the individual chunks
        big_chunk_size = chunk_size*((y_chunks)+1)
        dont_add = False
        if(direction==0):
            for line in data:
                if(line==[]):
                    if not dont_add:
                        new_data.append([])
                    continue
                reflected_line = line.copy()
                reflected_line[1] = float(reflected_line[1])
                if(reflected_line[1] == y_neg_bound):
                    dont_add = True
                    continue
                else:
                    dont_add = False
                    reflected_line[1] = -reflected_line[1]
                    reflected_line[1] = format(reflected_line[1],'.6f')
                    new_data.append(reflected_line)
            new_data.append([])
            for i in range(0,num_chunks):
                for j in range(0,chunk_size):
                    data.insert((int(i/y_chunks))*(chunk_size*y_chunks)+int((i/y_chunks))*(chunk_size*(y_chunks+1))+j,new_data[(i*chunk_size)+j])
        elif(direction==1):
            for line in data:
                if(line==[]):
                    if not dont_add:
                        new_data.append([])
                    continue
                reflected_line = line.copy()
                reflected_line[1] = float(reflected_line[1])
                if(reflected_line[1] == y_pos_bound):
                    dont_add = True
                    continue
                else:
                    dont_add = False
                    reflected_line[1] = y_pos_bound + (y_pos_bound-reflected_line[1])
                    reflected_line[1] = format(reflected_line[1],'.6f')
                    new_data.append(reflected_line)
            data.append([])
            for i in range(0,num_chunks):
                for j in range(0,chunk_size):
                    data.insert((int(i/y_chunks))*(chunk_size*y_chunks)+int((i/y_chunks))*(chunk_size*(y_chunks+1))+j+big_chunk_size,new_data[(i*chunk_size)+j])
    elif(axis==2):
        num_chunks = int((((z_pos_bound-z_neg_bound)/increment))*(((y_pos_bound-y_neg_bound)/increment)+1)) #number of individual blocks
        z_chunks = int(((z_pos_bound-z_neg_bound)/increment))
        y_chunks = int(((y_pos_bound-y_neg_bound)/increment)+1)
        chunk_size = int(((x_pos_bound-x_neg_bound)/increment)+2) #size of the individual chunks
        dont_add = False
        if(direction==0):
            for line in data:
                if(line==[]):
                    if not dont_add:
                        new_data.append([])
                    continue
                reflected_line = line.copy()
                reflected_line[2] = float(reflected_line[2])
                if(reflected_line[2] == z_neg_bound):
                    dont_add = True
                    continue
                else:
                    dont_add = False
                    reflected_line[2] = -reflected_line[2]
                    reflected_line[2] = format(reflected_line[2],'.6f')
                    new_data.append(reflected_line)
            new_data.append([])
            for k in range(0,z_chunks):
                for i in range(0,int(y_chunks/2)):
                    swap_with = y_chunks-1-i
                    for j in range(0,chunk_size):
                        temp1 = new_data[(k*(chunk_size*y_chunks))+(chunk_size*i)+j].copy()
                        temp2 = new_data[(k*(chunk_size*y_chunks))+(chunk_size*swap_with)+j].copy()
                        new_data[(k*(chunk_size*y_chunks))+(chunk_size*i)+j] = temp2
                        new_data[(k*(chunk_size*y_chunks))+(chunk_size*swap_with)+j] = temp1
            print(new_data)
            for i in range(0,num_chunks):
                for j in range(0,chunk_size):
                    data.insert(j,new_data[j+(i*chunk_size)])
        elif(direction==1):
            for line in data:
                if(line==[]):
                    if not dont_add:
                        new_data.append([])
                    continue
                reflected_line = line.copy()
                reflected_line[2] = float(reflected_line[2])
                if(reflected_line[2] == z_pos_bound):
                    dont_add = True
                    continue
                else:
                    dont_add = False
                    reflected_line[2] = z_pos_bound + (z_pos_bound-reflected_line[2])
                    reflected_line[2] = format(reflected_line[2],'.6f')
                    new_data.append(reflected_line)
            data.append([])
            new_data.append([])
            for k in range(0,z_chunks):
                for i in range(0,int(y_chunks/2)):
                    swap_with = y_chunks-1-i
                    for j in range(0,chunk_size):
                        temp1 = new_data[(k*(chunk_size*y_chunks))+(chunk_size*i)+j].copy()
                        temp2 = new_data[(k*(chunk_size*y_chunks))+(chunk_size*swap_with)+j].copy()
                        new_data[(k*(chunk_size*y_chunks))+(chunk_size*i)+j] = temp2
                        new_data[(k*(chunk_size*y_chunks))+(chunk_size*swap_with)+j] = temp1
            for i in range(0,num_chunks):
                for j in range(0,chunk_size):
                    data.insert((chunk_size*(num_chunks+y_chunks))+j,new_data[j+(i*chunk_size)])
    return data

def write_to_file(file_name,data):
    str_to_append = ""
    for line in data:
        if(line==[]):
            str_to_append += "\n"
        else:
            line = "   ".join(line)
            str_to_append += line
            str_to_append += "\n"
    append_data = open(file_name,"w+")
    append_data.write(str_to_append)
    append_data.close()

while(True):
    data = storeData(file_name)
    last_el = data[len(data)-1]
    secondlast_el = data[len(data)-2]
    first_el = data[0]
    x_pos_bound = float(last_el[0])
    x_neg_bound = float(first_el[0])
    if(num_dims>1):
        y_neg_bound = float(first_el[1])
        y_pos_bound = float(last_el[1])
        if(num_dims>2):
            z_neg_bound = float(first_el[2])
            z_pos_bound = float(last_el[2])
    increment = x_pos_bound - float(secondlast_el[0])
    print(x_pos_bound)
    print(y_pos_bound)
    print(z_pos_bound)
    print(x_neg_bound)
    print(y_neg_bound)
    print(z_neg_bound)
    print(increment)
    axis_inp = input("Reflect on which axis? (x/y/z): ")
    direction_inp = input("Reflect in which direction? (+/-): ")
    reflect_axis = -1
    reflect_dir = -1
    if(axis_inp=="x"):
        if(direction_inp=="-"):
            reflect_axis = 0
            reflect_dir = 0
        elif(direction_inp=="+"):
            reflect_dir = 1
            reflect_axis = 0
        else:
            print("Bad direction input. Restarting.")
    elif(axis_inp=="y"):
        if(direction_inp=="-"):
            reflect_dir = 0
            reflect_axis = 1
        elif(direction_inp=="+"):
            reflect_dir = 1
            reflect_axis = 1
        else:
            print("Bad direction input. Restarting.")
    elif(axis_inp=="z"):
        if(direction_inp=="-"):
            reflect_dir = 0
            reflect_axis = 2
        elif(direction_inp=="+"):
            reflect_dir = 1
            reflect_axis = 2
        else:
            print("Bad direction input. Restarting.")
    else:
        print("Bad axis input. Restarting.")
        continue
    more_data = reflect(reflect_axis,reflect_dir,data)
    write_to_file(out_file_name,more_data)
    file_name = out_file_name
    go = input("Press y to continue, anything else to exit: ")
    if(go is not "y"):
        break

