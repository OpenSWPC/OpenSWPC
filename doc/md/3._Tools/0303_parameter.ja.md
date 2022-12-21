# パラメタ作成支援

## `fdmcond.x` 

差分法の空間・時間グリッド間隔は，主としてCFL安定条件に律速される．また，震源で輻射される地震動の特徴的周波数の上限は波長条件（地震波の波長が5$\sim$10グリッドより大きくなくてはならない）によって制限される．
前者が満たされないと計算は安定に実行されず発散し，後者の場合には発散こそしないものの数値分散にともなう誤差が大きくなる．

fdmcond.x
はこれらのパラメタを対話的に求めるツールである．空間・時間のグリッド間隔，最大周波数fmax，震源におけるライズタイムTr，地震波速度の最小最大値vmin,
vmaxなどの幾つかの組み合わせを与えると，残りのパラメタの条件を提示してくれる．

``` 
$  ./fdmcond.x 
 
----------------------------------------------------------------------
                           FDM CONDITION                           
----------------------------------------------------------------------
 
 
  Model Dimension ? --> 3
   2) 2D
   3) 3D
 
 
 Source Type ? --> 3
   1) Triangle
   2) Herrmann
   3) Kupper
 
 
 Parameter Combination ? --> 5
   1) dh   (space grid),  fmax (max freq.),  vmax (max vel.)
   2) dh   (space grid),  Tr   (rise time),  vmax (max vel.)
   3) dh   (space grid),  fmax (max freq.),  dt (time grid)
   4) dh   (space grid),  Tr   (rise time),  dt (time grid)
   5) dh   (space grid),  vmin (min vel.),   vmax (max vel.)
   6) dh   (space grid),  vmin (min vel.),   dt (time grid) 
   7) fmax (max freq.) ,  vmax (max vel.),   dt (time grid)
   8) Tr   (rise time) ,  vmax (max vel.),   dt (time grid)
   9) vmin (min vel.)  ,  vmax (max vel.),   dt (time grid)
 
 
 Assumed Parameters: 
   dx     =   0.25
   dy     =   0.25
   dz     =   0.25
   vmin   =   0.3
   vmax   =   8.0
 
 Derivaed Parameters: 
   dt    <=   0.01546
   fmax  <=   0.17143
   Tr    >=  13.41667
 
 
```

## `mapregion.x`

計算領域は`clon`, `clat`, `phi`, `xbeg`, `ybeg`, `nx`, `ny`, `dx`,
`dy`などのパラメタによって自動的に決定される．
`mapregion.x`は，の入力パラメタファイルを読みこみ，対応する地図領域範囲を出力するプログラムである．`xbeg`, `ybeg`, `xend`, `yend`で定義される領域の外側線の緯度経度位置を出力する．

``` bash
$ mapregion.x -i input.inf -o region.dat
```

`-o`オプションが与えられない場合には標準出力に出力を行う．
また，本コードは$N_M=3$を仮定した場合の予想メモリ利用量をターミナルに出力する．

## `mapregion.gmt`

`mapregion.x`を呼び出して地図領域を描画する`GMT`スクリプト．GMTバージョン5以降が必要である．
