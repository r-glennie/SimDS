# Python script reads MCDS output files and returns csv file 
# with parameter estimates for each simulated output file. 

# user options to set: 
results_filename = "results.txt" 
no_of_files = 100 
mcds_output_directory = "./" 

# csv file that estimates are to be written to 
results = open(results_filename,'w')
adjust=0 

# headers 
results.write("A1, A2, A3, density, se \n")
for i in range(1,no_of_files+1):
    # assumed that output files are named "output_<simulation_number>.txt" 
    filename = mcds_output_directory + "output_"+str(i)+".txt"
    outfile = open(filename,'r')
    for line in outfile:
        next=0
        if(len(line)>=4):
            if(line[3]==" " and line[4]=="D" and line[5]==" "):
                results.write(line[14:20]+",")
                results.write(line[27:34]+"\n")
            if (line[4:9]=="A( 1)"):
                results.write(line[15:20]+",")
            if (line[4:9]=="A( 2)"):
                results.write(line[15:20]+",")
                next=1
            if (adjust==1):
                if (line[4:9]=="A( 3)"):
                    results.write(line[13:20]+",")    
                else:
                    results.write("0 ,")
                adjust=0
            if (next==1): 
                adjust=1

    outfile.close()
results.close()
