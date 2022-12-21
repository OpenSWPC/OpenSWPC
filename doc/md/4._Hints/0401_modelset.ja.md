# 計算モデル設定

特に3次元計算では，計算規模はメモリサイズによって制限される．本コードが必要とするメモリ量の概算値は，混合精度かつ`NM=3`の場合，

\begin{align}
\begin{split}
    &m_\text{MP} = 116 + 24 \mathtt{NM} = 188 \quad (\mathtt{NM}=3) \quad \text{byte} \\
\end{split}\end{align}

である．ただし，これは3次元以上の配列の容量の総和を取ったもので，かつ吸収境界条件部分は考慮していない．
これに$n_x \times n_y \times n_z$を乗ずれば必要とする合計メモリ量の概算値が得られる．
あとは各計算ノードでの要求メモリ量が利用可能な資源量を超えないように設定すれば良い．

計算時間は，1秒あたりに1CPUコアが処理できる時間・空間グリッド数（$\mathtt{n_G}$）から推定できる．次表にいくつかのマシンでの推定結果を示す．全体の計算時間は，

\begin{align}
    t_\text{comp} = \frac{ {\tt nx}\times {\tt ny}\times {\tt nz}}{ \mathtt{n_G} \times \tt nproc} \times {\tt nt}
    \quad \text{[s]}\end{align}

により概算できる．計算時間がシステムの制限値を超えるようであれば，モデルサイズを変更するか，あるいはチェックポイントリスタート機能を利用する必要がある．


| マシン名         | CPU                           | Core数 | $\mathtt{n_G}$      |
| ---------------- | ----------------------------- | ------ | ------------------- |
| Mac Pro 2010     | Intel Xeon X5670 2.93GHz      | 6      | $6.7 \times 10^{6}$ |
| EIC2015          | Intel Xeon E5-2680 v3 2.5 GHz | 12     | $7.0 \times 10^{6}$ |
| 地球シミュレータ | NEC SX-ACE                    | 4      | $57 \times 10^{6}$  |