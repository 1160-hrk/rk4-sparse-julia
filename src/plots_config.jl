#!/usr/bin/env julia
"""
シンプルなPlots.jl設定ライブラリ
論文・プレゼン用のグラフ設定を簡単に変更
"""

module PlotsConfig

using Plots

# プリセット設定
const PRESETS = Dict(
    "paper" => Dict(
        :titlefontsize => 12,
        :guidefontsize => 11,
        :tickfontsize => 10,
        :legendfontsize => 10,
        :linewidth => 1.5,
        :size => (960, 540),  # 16:9ピクセル
        :dpi => 150
    ),
    "presentation" => Dict(
        :titlefontsize => 18,
        :guidefontsize => 16,
        :tickfontsize => 14,
        :legendfontsize => 12,
        :linewidth => 3.0,
        :size => (960, 540),
        :dpi => 100
    ),
    "presentation_large" => Dict(
        :titlefontsize => 28,
        :guidefontsize => 24,
        :tickfontsize => 20,
        :legendfontsize => 18,
        :linewidth => 4.0,
        :size => (960, 540),
        :dpi => 100
    )
)

"""
apply_style(preset::String; kwargs...)

スタイルプリセットを適用

Parameters
----------
preset : "paper", "presentation", "presentation_large"のいずれか
kwargs : 追加カスタマイズ設定
"""
function apply_style(preset::String = "presentation"; kwargs...)
    if !haskey(PRESETS, preset)
        error("不明なプリセット: $preset. 利用可能: $(keys(PRESETS))")
    end
    settings = deepcopy(PRESETS[preset])
    merge!(settings, Dict(kwargs))
    Plots.default(; settings...)
    _apply_common_settings()
end

"""
reset()

デフォルト設定に戻す
"""
function reset()
    Plots.default()
end

"""
set_figsize(width::Real, height::Real)

図のサイズを指定
"""
function set_figsize(width::Real, height::Real)
    Plots.default(size = (round(Int, width), round(Int, height)))
end

"""
list_presets()

利用可能なプリセット一覧を返す
"""
list_presets() = collect(keys(PRESETS))

"""
with_temp_style(preset::String; kwargs...)(f)

一時的にスタイルを適用するdoブロック

Example
-------
with_temp_style("paper") do
    plot(rand(10), rand(10))
end
"""
function with_temp_style(preset::String; kwargs...)
    orig = deepcopy(Plots.plotattr())
    return function(f::Function)
        try
            apply_style(preset; kwargs...)
            f()
        finally
            Plots.default(; orig...)
        end
    end
end

function _apply_common_settings()
    Plots.default(
        gridalpha = 0.3,
        # background_color_inside = :transparent,
        # background_color_outside = :transparent,
        # legend_background_color = :transparent,
        fg_color_subplot = :transparent,
        guidefontcolor = :black,
        tickfontcolor = :black,
        framestyle = :box,
        margin = 10Plots.mm,  # 全方向に10mmのマージンを追加
        bottom_margin = 12Plots.mm,  # 下部に追加マージン（x軸ラベル用）
        left_margin = 12Plots.mm     # 左側に追加マージン（y軸ラベル用）
    )
end

# モジュール読み込み時にpresentationスタイルを適用
apply_style("presentation")

# 使用例
if abspath(PROGRAM_FILE) == @__FILE__
    println("利用可能なプリセット:")
    for p in list_presets()
        println("  - $p")
    end

    println("\n使用例:")
    println("""  using PlotsConfig""")
    println("""  PlotsConfig.apply_style("paper")  # 論文用""")
    println("""  PlotsConfig.set_figsize(1200, 675)  # サイズ変更""")
    println("""  PlotsConfig.with_temp_style("paper") do""")
    println("""      plot(rand(10), rand(10))  # 一時的に適用""")
    println("""  end""")
end

end # module
