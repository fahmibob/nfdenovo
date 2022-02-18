nextflow.enable.dsl = 2
params.pair= "$baseDir/input_data/data1_{fwd,rev}.txt"
params.remove= true

workflow {
  channel.fromFilePairs(params.pair, flat: true)
    .map {data -> tuple(data[1].parent.getName(), data[0], data[1], data[2] )}
    .set{pair_ch}

  pair_ch.view()
  //test.out.toList().view()
    //pair_ch.view()
}
Workflow pair {
  
}
/*
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
*?
