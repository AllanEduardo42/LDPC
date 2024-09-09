using Plots
using LaTeXStrings

plotlyjs()
plot(
    1:4,
    [[1,4,9,16]*10000, [0.5, 2, 4.5, 8]],
    labels = [L"\alpha_{1c} = 352 \pm 11 \text{ km s}^{-1}";
              L"\beta_{1c} = 25 \pm 11 \text{ km s}^{-1}"] |> permutedims,
    xlabel = L"\sqrt{(n_\text{c}(t|{T_\text{early}}))}",
    ylabel = L"d, r \text{ (solar radius)}",
    yformatter = :plain,
    extra_plot_kwargs = KW(
        :include_mathjax => "cdn",
        :yaxis => KW(:automargin => true),
        :xaxis => KW(:domain => "auto")
   ),
)