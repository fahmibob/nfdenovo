nextflow.enable.dsl=2

params.jobs="single"
params.threads=5
params.single=false
params.pair=false

params.fasta=["fasta","fa","fas","fna"]

params.kmer=35
params.avg_ins=200
params.reverse_seq=0
params.asm_flags=3
params.minContigLen=150
params.soap31="SOAPdenovo-Trans-31mer"
params.soap127="SOAPdenovo-Trans-127mer"
params.bbreformat="$baseDir/tools/bbmap/reformat.sh"


workflow {
  if (params.single) {
    channel.fromPath(params.single)
      .filter( ~/.*\d+.fq.gz/ )
      .map {file -> tuple(file.parent.getName(), file.simpleName, file)}
      .set{single_ch}

    single_assembly(single_ch)
  } else if (params.pair) {
    channel.fromFilePairs(params.pair, flat: true)
      .map {data -> tuple(data[1].parent.getName(), data[0], data[1], data[2] )}
      .set{pair_ch}

    pair_assembly(pair_ch)
  }
}

process single_assembly {
  publishDir "$baseDir/output/${params.jobs}_${params.kmer}/${speciesName}", mode: 'copy'
  errorStrategy 'ignore'
  maxForks 4

  input:
  tuple val(speciesName), val(sampleName), path(single)

  output:
  path("*.fasta")

  script:
  if (single.getExtension() == "gz") {
    if (single.getBaseName().lastIndexOf('.').with {it != -1 ? single.getBaseName()[it+1..-1] : single.getBaseName()} in params.fasta) {
      configDataType="f"
      data_type="fasta"
    } else {
      configDataType="q"
      data_type="fastq"
    }
    single_data=single.getBaseName()
  } else {
    if (single.getExtension() in params.fasta) {
      configDataType="f"
      data_type="fasta"
    } else {
      configDataType="q"
      data_type="fastq"
    }
    single_data=single
  }
  """

  if [[ $single == *.gz ]]
  then
    pigz --fast -df $single
  fi
  echo "[LIB]" >> $sampleName'.config'
  echo "avg_ins="$params.avg_ins >> $sampleName'.config'
  echo "reverse_seq="$params.reverse_seq >> $sampleName'.config'
  echo "asm_flags="$params.asm_flags >> $sampleName'.config'
  echo $configDataType"="$single_data >> $sampleName'.config'

  if [[$params.kmer < 32 ]]
  then
    $params.soap31 all -s $sampleName'.config' -L $params.minContigLen \
    -K $params.kmer -o $sampleName'_'$params.kmer
  else
    $params.soap127 all -s $sampleName'.config' -L $params.minContigLen \
    -K $params.kmer -o $sampleName'_'$params.kmer
  fi
  $params.bbreformat in=$sampleName'_'$params.kmer'.scafSeq' out=$sampleName'_'$params.kmer'.fasta' minlength=200

  rm $sampleName'_'$params.kmer'.scafSeq'
  rm -rf $single_data
  """
}

process pair_assembly {
  publishDir "$baseDir/output/${params.jobs}_${params.kmer}/${speciesName}", mode: 'copy'
  errorStrategy 'ignore'
  maxForks 5

  input:
  tuple val(speciesName), val(sampleName), path(forward), path(reverse)

  output:
  //tuple val(speciesName), val(sampleName), path(data_forward), path(data_reverse)
  //tuple path("sampled*_1*"), path("sampled*_2*")
  //tuple path("*.contig"),
  path("*.fasta")

  script:
  if (forward.getExtension() == "gz") {
    if (forward.getBaseName().lastIndexOf('.').with {it != -1 ? forward.getBaseName()[it+1..-1] : forward.getBaseName()} in params.fasta) {
      configDataType_fwd="f1"
      configDataType_rev="f2"
      data_type="fasta"
    } else {
      configDataType_fwd="q1"
      configDataType_rev="q2"
      data_type="fastq"
    }
    data_forward=forward.getBaseName()
    data_reverse=reverse.getBaseName()
  } else {
    if (forward.getExtension() in params.fasta) {
      configDataType_fwd="f1"
      configDataType_rev="f2"
      data_type="fasta"
    } else {
      configDataType_fwd="q1"
      configDataType_rev="q2"
      data_type="fastq"
    }
    data_forward=forward
    data_reverse=reverse
  }
  """
  if [[ $forward == *.gz ]]
  then
    pigz --fast -df $forward $reverse
  fi
  echo "[LIB]" >> $sampleName'.config'
  echo "avg_ins="$params.avg_ins >> $sampleName'.config'
  echo "reverse_seq="$params.reverse_seq >> $sampleName'.config'
  echo "asm_flags="$params.asm_flags >> $sampleName'.config'
  echo $configDataType_fwd"="$data_forward >> $sampleName'.config'
  echo $configDataType_rev"="$data_reverse >> $sampleName'.config'

  if [[ $params.kmer < 32 ]]
  then
    $params.soap31 pregraph -s $sampleName'.config' \
    -K $params.kmer -p $params.threads -o initial_$sampleName
    $params.soap31 contig -g initial_$sampleName
  else
    $params.soap127 pregraph -s $sampleName'.config' \
    -K $params.kmer -p $params.threads -o initial_$sampleName
    $params.soap127 contig -g initial_$sampleName
  fi

  cat $data_forward | head -n 800000 > sampled_$sampleName'_1.'$data_type
  cat $data_reverse | head -n 800000 > sampled_$sampleName'_2.'$data_type

  bwa index initial_$sampleName'.contig'
  bwa mem initial_$sampleName'.contig' sampled_$sampleName'_1.'$data_type sampled_$sampleName'_2.'$data_type | samtools view -Sb > $sampleName'.bam'
  samtools view $sampleName'.bam' | cut -f9 > inserts_$sampleName
  ruby -r "$baseDir/bin/module.rb" -e "avgInsertBam('inserts_$sampleName')" > avginsert_$sampleName
  sed -i -e 's/$params.avg_ins/'"\$(grep '[[:digit:]]' avginsert_$sampleName)/g" $sampleName'.config'

  if [[$params.kmer < 32 ]]
  then
    $params.soap31 all -s $sampleName'.config' -L $params.minContigLen \
    -K $params.kmer -p $params.threads -o $sampleName'_'$params.kmer
  else
    $params.soap127 all -s $sampleName'.config' -L $params.minContigLen \
    -K $params.kmer -p $params.threads -o $sampleName'_'$params.kmer
  fi
  $params.bbreformat in=$sampleName'_'$params.kmer'.scafSeq' out=$sampleName'_'$params.kmer'.fasta' minlength=200

  rm -rf initial_$sampleName*
  rm -rf $sampleName'.bam'
  rm -rf sampled_$sampleName'_1.'$data_type sampled_$sampleName'_2.'$data_type
  rm $sampleName'_'$params.kmer'.scafSeq'
  rm -rf $data_forward $data_reverse
  """
}
