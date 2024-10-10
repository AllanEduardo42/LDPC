# a = collect(1.0:1500.0)

# for i=1:10
#     idx = rand(1:1500)
#     a[idx] = randn()
# end

a = randn(1500)

b = copy(a)
@time sort!(b);