# 相反定理モード

観測点から点震源による地震波を励起し，複数の震源位置を仮想観測点とした波形を計算する．
これらは，相反定理によって，各震源位置から観測点までの計算波形に相当する．これらは，震源時間関数が充分に短ければGreen関数として扱うことができる．

本モードで計算されるのは，震源位置$\mathbf{\xi}$から観測点位置$\mathbf{r}$へのGreenテンソルを$G_{ij}\left( \mathbf{r}, t; \mathbf{\xi}\right)$に，計算時に仮定する震源時間関数を$s(t)$がたたみこまれた

\begin{align}
\begin{split}
  G^{M1}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \frac{\partial G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_x} \ast s(t) =   \frac{\partial G_{ix} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_x} \ast s(t) 
  \\
  G^{M2}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \frac{\partial G_{iy} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_y} \ast s(t) 
   =     \frac{\partial G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_y} \ast s(t)    
  \\
  G^{M3}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \frac{\partial G_{iz} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_z} \ast s(t) 
  =      \frac{\partial G_{iz} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_z} \ast s(t) 
  \\
  G^{M4}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \left( 
    \frac{\partial G_{iy} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_z} + 
    \frac{\partial G_{iz} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_y} \right) \ast s(t) 
    \\
  &= \left( 
    \frac{\partial G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_z} + 
    \frac{\partial G_{iz} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_y} \right) \ast s(t) 
  \\
  G^{M5}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \left(
    \frac{\partial G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_z} + 
    \frac{\partial G_{iz} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_x}\right) \ast s(t)  
    \\& = \left(
    \frac{\partial G_{ix} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_z} + 
    \frac{\partial G_{iz} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_x}\right) \ast s(t) 
  \\
  G^{M6}_{i} \left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv \left(
    \frac{\partial G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_y} + 
    \frac{\partial G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_x} \right) \ast s(t) 
    \\
  &= \left(
    \frac{\partial G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) }{\partial \xi_y} + 
    \frac{\partial G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) }{\partial \xi_x} \right) \ast s(t)
 \end{split} 
\end{align}

である．また，オプションで実体波Green計算

\begin{align}
\begin{split}
G^{B1}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv  G_{ix} \left( \mathbf{r}, t; \mathbf{\xi}\right) \ast s(t)
  =       G_{ix} \left( \mathbf{\xi}, t; \mathbf{r}\right) \ast s(t)  
  \\
G^{B2}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv G_{iy} \left( \mathbf{r}, t; \mathbf{\xi}\right) \ast s(t) 
  =      G_{iy} \left( \mathbf{\xi}, t; \mathbf{r}\right) \ast s(t) 
\\
  G^{B3}_{i}\left( \mathbf{r}, t; \mathbf{\xi}\right) 
  &\equiv G_{iz} \left( \mathbf{r}, t; \mathbf{\xi}\right)\ast s(t)
  =      G_{iz} \left( \mathbf{\xi}, t; \mathbf{r}\right)\ast s(t)
\end{split}
\end{align}

も計算される．

パラメタファイルで指定された観測点 `green_stnm` から，`green_cmp` 成分の実体力による地震波を，時間長さ`green_trise` の震源時間関数で励起する．
そのため，全グリーンテンソルを求めるには，`green_cmp=’x’, ’y’, ’z’` の3種類の計算を行う必要がある．

震源要素はデカルト座標系もしくは緯度経度および深さで指定し，それぞれの震源要素に対して整数の `gid` を付与する．
結果は，`(odir)/green/(gid)` ディレクトリに，`(title)__(green_cmp)__mij__.sac` 
あるいは`(title)__(green_cmp)__fi__.sac` というファイル名のSACファイルとして格納される．

`gid` は必ずしも連番である必要はない．SACの規約波形と直接比較するため，出力波形には $10^9$ が乗じられている．また，$z$は下が正になる座標系を取っているが，Green関数として出力される波形は，多くの場合の観測に準じ，上が正になるように変換されている．しかし，モーメントテンソル応答の計算の深さ方向微分は計算座標系の下が正の定義のまま行われている．

!!! Info "Parameters"
    `green_mode`
    :Green関数モードのON/OFFを設定する．`.true.`でONになり，通常の震源に関する指定は無視される．

    `green_stnm`
    : 受信観測点名．この観測点が別途与える観測点リストに含まれている必要がある．

    `green_cmp`
    : Green関数を計算する成分．受信点における観測成分に相当する．`’x’`, `’y’`, `’z’`のいずれかを指定．

    `green_trise`
    : Green関数にたたみ込まれる震源時間関数のライズタイム

    `green_bforce`
    : `.true.`のとき実体波のGreen関数も計算する．デフォルトは`.false.`

    `green_fmt`
    : Green関数の震源要素位置ファイルリストのフォーマット．`’xyz’`（デカルト座標；デフォルト）か `’llz’`（緯度・経度・深さ）のどちらかを指定する．

    `green_maxdist`
    : 震源要素位置と観測点の間の震央距離がこの指定パラメタ（km単位）以内の場合にのみGreen関数を計算する．デフォルトでは常に計算する．

    `fn_glst`
    : Green関数震源要素位置ファイル名．フォーマットは表の通り．

    `stftype`
    : 震源時間関数種別．指定法はモーメントテンソル震源に同じ．

    `ntdec_w`
    : Green関数波形出力の時間間引き間隔．通常計算の波形出力時の指定に同じ．

震源要素位置ファイルは`green_fmt`の値に応じて以下の書式を取る．

 | `green_fmt` | 書式 |     |     |       |
 | ----------- | ---- | --- | --- | ----- |
 | `’xyz’`     | x    | y   | z   | `gid` |
 | `’llz’`     | lon  | lat | z   | `gid` |
 
