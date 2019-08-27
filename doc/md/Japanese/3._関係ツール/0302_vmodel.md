# 構造モデルの確認

## `qmodel_tau.x`

`OpenSWPC`の入力パラメタファイルを元に，$Q^{-1}$値と実体波の分散性の周波数依存性を計算する．

``` bash
$ qmodel_tau.x -nm [nm] -i [prm_file] -f0 [min_freq] -f1 [max_freq] -nf [ngrid]
```

周波数`min_freq`から`max_freq`までの区間を対数等間隔に`ngrid`に分割し，$Q^{-1}(f)$と実体波速度の分散性を出力する．地震波速度は基準周波数で1に規格化されている．粘弾性体に関わるパラメタはの入力`prm_file`から読み込むが，コードの埋め込みパラメタである`NM`は別途オプションで与える必要がある．

## `grdsnp.x`

構造モデルgrdファイルから，`OpenSWPC`の入力パラメタに基づいた座標変換により，計算座標系の(x, y, depth)をアスキーで標準出力に出力する．計算モデル領域の切り出しの確認や，それをGMT等で可視化する際に用いる．

``` 
$ grdsnp.x -i [prm_file] -g [grd_file]
```
