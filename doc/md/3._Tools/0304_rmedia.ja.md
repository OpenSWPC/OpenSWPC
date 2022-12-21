# ランダム媒質生成

## `gen_rmed3d.x` 

3次元ランダム媒質を作成する．

```
gen_rmed3d.x  [-o outfile] [-nx nx] [-ny ny] [-nz nz] [-epsil epsilon]  [-kappa kappa] 
               [-dx dx] [-dy dy] [-dz dz] [-ax ax] [-ay ay] [-az az]  
               [-ptype ptype] {-seed seed_number}
```

`-o outfile`

:    
    出力ファイル名

`-nx nx`, `-ny ny` `-nz nz`

:    
    グリッドサイズ．2のべき乗であること．

`-epsil epsilon`

:    
    ランダム媒質のゆらぎRMS $\epsilon$

`-ax ax`, `-ay ay`, `-az az`

:    
    ランダム媒質の各方向の特徴的スケール．

`-dx dx -dy dy -dz dz`

:    
    空間グリッド間隔．シミュレーションで用いる間隔と同一でなければならない．

`-ptype ptype`

:    
    PSDFの種別．1:Gaussian, 2:Exponential, 3:von Kármán

`-kappa kappa`

:    
    von Kármán型のパラメータ$\kappa$

`-seed seed_number`

:    
    （オプション）乱数のシードを固定．整数値で与える．省略された場合は日付時刻に基づき自動生成される．

ランダム媒質は3次元の`NetCDF`フォーマットで作成される．一般的な可視化ソフト（たとえば
[`ParaView`](https://www.paraview.org)や`NetCDF`専用の可視化ソフト（たとえば[`Panoply`](https://www.giss.nasa.gov/tools/panoply/)）などで可視化したり内容を確認することができる．

!!! Quote "ランダム媒質の例"
    ![](../fig/rmedia_sample.png)
    `ParaView`による3次元ランダム媒質の可視化の例．

## `gen_rmed2d.x`

2次元ランダム媒質を作成する．使い方は`gen_rmed3d.x`とほぼ同じであり，y軸方向の入力が無視される．

