using DelimitedFiles

str = "C-RBP 0.85"

directory = "./Saved Data/join/"
directory2 = "./Saved Data/Artigo EbN0 576 3.5dB/"
directory3 = "./Saved Data/Artigo EbN0 576/"

FB = ["F","B"]

for j = 1:2
    for k in axes(modes_markers,1)
        x = readdlm(directory*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
        y = readdlm(directory2*FB[j]*"ER_"*str*".txt",'\t',Float64,'\n')
        z = [x y]
        open(directory3*FB[j]*"ER_"*str*".txt","w") do io
            writedlm(io,z)
        end
    end
end