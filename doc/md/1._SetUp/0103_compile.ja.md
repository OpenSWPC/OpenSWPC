# コンパイル

## make

`src/swpc_3d`, `src/swpc_psv`, `src/swpc_sh`, `src/tools`の下に`makefile`がある．
`src/shared`には各プログラムから共通に呼び出されるモジュールコードが格納されている．
各ソースコードディレクトリ（`src/swpc_*`）において`make`することで`bin`以下に実行バイナリが生成される．

## makefileの指定変数

makefileでは以下の変数を指定する必要がある．

| 変数     | 説明                                |
| -------- | ----------------------------------- |
| `FC`     | コンパイラ                          |
| `FFLAGS` | コンパイルオプション                |
| `NCFLAG` | `NetCDF` 利用スイッチ               |
| `NCLIB`  | `NetCDF` ライブラリ格納ディレクトリ |
| `NCINC`  | `NetCDF` includeディレクトリ        |
| `NETCDF` | `NetCDF` ライブラリのリンク         |


`NCFLAG = -D_NETCDF` が指定されることにより，`NetCDF`モジュールを利用したコンパイルが行われる．


さまざまな計算機環境でコンパイルするため，以下の`arch`オプションによりオプションを分岐させている: 

| archオプション名 | 対象                                   | NetCDFライブラリの場所  |
| ---------------- | -------------------------------------- | ----------------------- |
| mac-intel        | Mac OSX + Intel Compiler + OpenMPI     | `${HOME}/local`         |
| mac-gfortran     | Mac OSX + gfortran + Open MPI          | `/usr/local`            |
| eic              | EIC2015環境下のIntel Compiler          | `${HOME}/local`         |
| fx  (obsolete) | 富士通FX10, FX100 および京コンピュータ | `${HOME}/xlocal`        |
| es3 (obsolete)             | 地球シミュレータ                       | システム提供            |
| ubuntu-gfortran  | Ubuntu 16.04LTS + gfortran + Open MPI  | `apt`によるインストール |
| ofp (or oak)  (obsolete)| 東京大学情報基盤センター Oakforest-PACS | `module` コマンドを通じて自動指定 |
| obcx (obsolete) | 東京大学情報基盤センター Oakbridge-CX | `module` コマンドを通じて自動指定 |
| bdec-o |  東京大学情報基盤センター Wisteria-Odyssey/BDEC01 |  `module` コマンドを通じて自動指定 |
| mac-m1  | macOS + gfortran (Apple Silicon M1/M2 + Homebrew) | `/opt/homebrew/` |


たとえば`mac-gfortran`に相当する環境では，

```make
make arch=mac-gfortran
```

とすることで，その環境に適したコンパイルオプションが自動的に選択される．また，幾つかの環境では，

```make
make arch=eic debug=true
```

のように`debug=true`オプションを付けると，コンパイルオプション`FFLAGS`がデバッグに適したものに変更されるようになっている．これらの変数は
`src/shared/makefile.arch`と`src/shared/makefile-tools.arch`に定義されている．新たな環境を追加するには，これらのファイルにオプションを追記するのが簡単であろう．


## `NetCDF`の利用

本コードの速度構造入力の一部には`NetCDF`形式を採用しており，ライブラリとモジュール情報ファイルが必要である．
具体的には，

- `libnetcdf.*`:   `NetCDF`ライブラリファイル
- `libnetcdff.*`:  `NetCDF` Fortran用ライブラリファイル（`NetCDF` ver.4以降）
- `netcdf.mod`:    Fortranモジュール情報ファイル

がコンパイル時に必要となる．
ライブラリの拡張子は`*.a`（static）の場合と`*.so`（dynamic）の場合がある．
また，Fortran90の仕様により，`NetCDF.mod`はのコンパイルと同じコンパイラにより作成されていなければならない．
特にLinux環境等において，`yum, apt, brew`等のパッケージ管理システムにより導入した`NetCDF`では，`gfortran`以外のコンパイラを用いることができない．そのような場合は，別の場所に用の`NetCDF`を自力でコンパイルする必要がある．

## 埋め込みパラメータの調整

`OpenSWPC`の挙動は原則としてパラメータファイルで制御されるが，計算速度の向上のため，幾つかの変更頻度の低いと期待されるパラメタがソースコードに埋め込まれている．これらはいずれも`m_global.F90`で定義されている．以下の値を修正した場合には，コード全体の再コンパイル（`make clean; make`）が必要である．

!!! Info "Parameters"
    **`UC`**
    : 計算に用いる混合座標系（距離km, 密度g/cm$^3$, 地震波速度km/sなどの慣用単位）での結果をSI単位に換算する定数． 計算の単位系を変更する場合には調整が必要となる．デフォルト値は `1e-15` (3D) もしくは `1e-12` (2D)
    
    **`MP`**
    : 混合精度差分法の精度．`DP`の場合，必要な変数については倍精度で，それ以外の変数を単精度でそれぞれ計算する．これをSPに変更すると，計算すべてを単精度で実行するようになる．単精度計算では特に震源周辺で桁あふれに伴う不安定が起こる場合があるが，必要となるメモリ量は約2/3に小さくなり，計算も速い．デフォルト値は倍精度を意味するパラメタ`DP`．具体的な値はFortranコンパイラによって自動的に決定されるが，多くのコンパイラで`8`を取る．
    
    **`NM`**
    : メモリ変数を用いて計算する一般化Zener粘弾性体の並列数．1より大きな値の場合には，パラメタ`fq_*`に応じて指定された周波数範囲において周波数一定の$Q$値を取るように最適化される．一方，0の場合には完全弾性体の問題となり，内部減衰は考慮されない．特に3次元計算においては，計算に必要な時間やメモリ量がこのパラメタによって大きく変動する．