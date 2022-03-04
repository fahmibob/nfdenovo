
nextflow.enable.dsl = 2
params.pair= "$baseDir/input_data/data1_{fwd,rev}.txt"
params.data= "$baseDir/data/*{fwd,rev}.txt"
params.remove= true
params.folderToFilter = "$baseDir/filterfolder"
params.listToFilter = "$baseDir/filterdata.txt"

workflow {
  //pair data
  Channel.fromFilePairs(params.data, flat: true)
    .set{pair_ch}

  //filterFolder=Channel.fromPath(params.folderToFilter).toList()


  if (params.listToFilter){
    def filterList = new File(params.listToFilter).collect {it}
    pair_ch=pair_ch.map { if (!(it[0].toString() in filterList)){it}}
  }
/*
  if (params.folderToFilter){
    def filterFolder = new FileNameByRegexFinder().getFileNames(params.folderToFilter, /.*\.txt/)
    println filterFolder
    filterFolder=filterFolder.collect {it.substring(it.lastIndexOf("/")+1, it.indexOf("."))}
    pair_ch=pair_ch.map { if (!(it[0].toString() in filterFolder)){it}}
  }
*/



/*
  if (params.folderToFilter){
    pair_ch=pair_ch.map { if (!(it[0].toString() in filterList)){
        it[0]
      } }
  }
*/

  //filterFolder=Channel.fromPath(params.folderToFilter)
  //listToFilter

  //folderToFilter

  //filterList.view()
  pair_ch.view()
}

/*

Channel.fromPath(params.data)
  .map {data -> data.getFileName()}
  .set{test_data}

retrieved= Channel.fromPath(params.filterFolder)
  .map {data -> tuple(data.getFileName(),data.getBaseName())}

test_data2=test_data.join(retrieved, remainder: true)

test_data2=test_data2.filter{it[1] == null}.map {data -> data[0]}

//data22=file(params.filterdata, checkIfExists: true)
ch_filterdata= Channel.from(file(params.filterdata).readLines())
  .map { row -> file(row, checkIfExists: true) }
  .map {data -> tuple(data.getFileName(),data.getBaseName())}


test_data3=test_data2.join(ch_filterdata, remainder: true)
test_data3.view()

channel.fromFilePairs(params.pair, flat: true)
  .map {data -> tuple(data[1].parent.getName(), data[0], data[1], data[2] )}
  .set{pair_ch}

pair_ch.view()
process test {
  input:
  tuple val(x), val(y), path(fwd), path(rev)

  output:
  tuple val(x), val(y), path(fwd), path(rev)
  """
  """
}

if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}
  ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()
*/
