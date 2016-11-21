results = open("results.txt",'w')

no_of_files = 100

adjust=0 

results.write("A1, A2, A3, density, se \n")
for i in range(1,no_of_files+1):
    filename = "./SavedOutputs/output"+str(i)+".txt"
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
