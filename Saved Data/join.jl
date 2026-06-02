using DelimitedFiles

directory1 = "./Saved Data/2026-05-08 22:11 5GNR 576 1|2 3.5dB 4.0dB 4.5dB 5.0dB 5.5dB/"
directory2 = "./Saved Data/2026-05-09 08:41 5GNR 576 1|2 6.0dB/"
directory3 = "./Saved Data/2026-05-08 22:11 5GNR 576 1|2 3.5dB 4.0dB 4.5dB 5.0dB 5.5dB 6.0dB/"
# mkdir(directory3)

FB = ["F","B"]
modes =["Flooding"                       
        "LBP"                  
        "RBP"                  
        "RD-RBP"                               
        "NW-RBP"                        
        "SVNF"                 
        "List-RBP"             
        "UBP-RBP"              
        "C-RBP"                
        "C&R-RBP"              
        "C&DR-RBP"]  
for j =1:2
    for mode in modes
        x = readdlm(directory1*FB[j]*"ER_"*mode*".txt",'\t',Float64,'\n')
        y = readdlm(directory2*FB[j]*"ER_"*mode*".txt",'\t',Float64,'\n')
        z = [x y]
        open(directory3*FB[j]*"ER_"*mode*".txt","w") do io
            writedlm(io,z)
        end
    end
end