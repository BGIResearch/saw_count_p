
# change the libPath to your specific path
libPath="/root/lib"

cmakePath="$libPath/cmake-3.17.2"
gccPath="$libPath/gcc-9.1.0"
CLI11Path="$libPath/CLI11-1.9.0"
spdlogPath="$libPath/spdlog-1.5.0"
libdeflatePath="$libPath/libdeflate-1.5"
htslibPath="$libPath/htslib-1.14"
yggPath="$libPath/ygg-master"
doctestPath="$libPath/doctest-2.3.7"
libxPath="$libPath/libx-1.1"
rwqueuePath="$libPath/readerwriterqueue-1.0.5"
hdf5Path="$libPath/hdf5-1.12.1"

binPath="/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/bin"
export PATH="$gccPath/bin:$cmakePath/bin:$binPath"

export LD_LIBRARY_PATH="$libdeflatePath/lib:$htslibPath/lib:$gccPath/lib64:$LD_LIBRARY_PATH"
export LIBRARY_PATH=$LD_LIBRARY_PATH

export C_INCLUDE_PATH="$hdf5Path/include:$rwqueuePath:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$libxPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$doctestPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$libdeflatePath/include:$yggPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$htslibPath/include:$gccPath/include:$C_INCLUDE_PATH"
export C_INCLUDE_PATH="$CLI11Path/include:$spdlogPath/include:$C_INCLUDE_PATH"
export CPLUS_INCLUDE_PATH=$C_INCLUDE_PATH

export CC="$gccPath/bin/gcc"
export CXX="$gccPath/bin/g++"

absPath(){
    relativePath=$1
    mkdir -p $relativePath
    if [ ${relativePath:0:1} == "/" ]
    then
        echo $relativePath
        return 0
    fi
    ap="$(cd $relativePath; pwd)"
    echo "$ap"
}

install(){
    libFile="$1"
    tlib="$installPath/lib"
    if [ -d "$tlib" ]
    then
        if [ -e "$libFile" ]
        then
            echo "Installing: $libFile into $tlib"
            cp $libFile $tlib
        else
            echo "[WARN] File not exists: $libFile"
        fi
    fi
}

srcPath="$(cd $(dirname $(dirname $0)); pwd)"
buildPath="build"
buildPath="$(absPath $buildPath)"
mkdir -p $buildPath
cd $buildPath
installPath="$buildPath/install"
mkdir -p $installPath/bin $installPath/lib

timeStart=$(date +%s)

test="OFF"
if [ $# == 1 ]
then
    if [ $1 == "ON" ]
    then
        test="ON"
    fi
fi

cmake $srcPath -DINSTALL_PATH=$installPath -DUNITTEST=$test -DHDF5PATH=$hdf5Path/lib

thread=$(grep -c ^processor /proc/cpuinfo)
make -j $thread #VERBOSE=1

if [[ $? == 0 ]];
then
    install $htslibPath/lib/libhts.so.3
    install $libdeflatePath/lib/libdeflate.so.0
    install $gccPath/lib64/libstdc++.so.6
    install $gccPath/lib64/libgomp.so.1
    install $gccPath/lib64/libgcc_s.so.1
    install $hdf5Path/lib/libhdf5.so.200
    mkdir -p $installPath/doc
    cp $srcPath/doc/Bam2Gem_workflow.html $installPath/doc

    cd $srcPath
    cp -R $installPath $srcPath
    echo -n "Success! "
else
    echo -n "Fail! "
fi

timeEnd=$(date +%s)
secs=$(($timeEnd - $timeStart))
mins=$(($secs/60))
secs=$(($secs%60))
echo "Cost: $mins Mins $secs Secs"
