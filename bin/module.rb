def percentile(values, percentile)
    values_sorted = values.sort
    k = (percentile*(values_sorted.length-1)+1).floor - 1
    f = (percentile*(values_sorted.length-1)+1).modulo(1)

    return values_sorted[k] + (f * (values_sorted[k+1] - values_sorted[k]))
end

def avgInsertBam(input)
  data=File.readlines(input)
  data_int=data.map{|x| x.gsub("\n","").to_i}
  data_rn=data_int.select{|x| x > 0.1}
  min=percentile(data_rn,0.15)
  max=percentile(data_rn,0.85)
  data_rn=data_rn.select{|x| x>=min && x <= max}
  sum_data_rn=data_rn.inject(0){|sum, x| sum + x}
  print sum_data_rn/data_rn.count
end

def assembly_stats(assemblyData)
  lengths=get_lengths(assemblyData).values
  filename=assemblyData.split("/")[-1]
  total = lengths.inject(:+)
  n2=total/2
  len = lengths.length
  sorted = lengths.sort
  sorted_rev = lengths.sort.reverse
  average = total / len
  median = len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2
  sum=0

  csum=sorted_rev.map {|x|sum += x}
  csum2=csum.select{|x| x >= n2}.min
  n50_idx=csum.index(csum2)
  stats=[filename, average, median, lengths, sorted_rev[n50_idx]]
  return stats
end

def get_lengths(input)
  inputData=File.read(input).split(">").drop(1)
  seqIdLen={}
  inputData.each do |indv|
    nl=indv.index("\n")
    header=indv[0..nl-1]
    header=header.split("\s")[0]
    seq=indv[nl+1..-1]
    seq.gsub!(/\n/,'')
    seqIdLen[header]=seq.length
  end
  return seqIdLen
end

def contig(assemblyData)
  ass=File.read(assemblyData).split(">").drop(1)
  contigId=[]
  ass.each do |indv|
    nl=indv.index("\n")
    header=indv[0..nl-1]
    contigId.push(header.split("\s")[0])
  end
  return contigId
end

def benchmark(blastOutData, assemblyData, idThres, lenThres, refData)
  parsedRef=get_lengths(refData)
  transcriptList=parsedRef.keys
  parsedassemblyID=contig(assemblyData)

  matchBlast=File.readlines(blastOutData).select{|x| x[0]!="#"}
  matchBlast.map!{|x| x.split("\t")}

  fileterd_id_len_matchBlast=[]
  filtered_id_matchBlast=matchBlast.select{|x| x[2].to_i >= idThres} #test
  filtered_id_matchBlast.select do |x|
      len=(x[9].to_i-x[8].to_i).abs().to_f

      if len/parsedRef[x[1]].to_f >= lenThres
        fileterd_id_len_matchBlast.push(x)
      end
  end
  truePositive=fileterd_id_len_matchBlast.map{|x| x[0]}
  truePositive=truePositive.uniq
  trueNegative=parsedassemblyID.-(truePositive)

  foundTranscript=fileterd_id_len_matchBlast.map{|x| x[1]}
  missingTranscript=transcriptList.-(foundTranscript.uniq)

  tp_count=truePositive.count
  tn_count=trueNegative.count
  ci=(tp_count.to_f/tn_count.to_f).round(2)
  accuracy=(tp_count.to_f/(tp_count+tn_count+missingTranscript.count).to_f).round(2)
  missingTranscriptratio=(missingTranscript.count.to_f/transcriptList.count.to_f).round(2)
  return ci, accuracy, missingTranscriptratio
end

def report_summary(assemblyData, blastOutData, refData, idThres, lenThres, coverage)
  stats=assembly_stats(assemblyData)
  benchmark_res=benchmark(blastOutData, assemblyData, idThres, lenThres, refData)

  puts "Information for assembly\t'#{stats[0]}'"
  puts "mean size\t#{stats[1]}"
  puts "median\t#{stats[2]}"
  puts "Largest\t#{stats[3].max}"
  puts "shortest\t#{stats[3].min}"
  puts "N50\t#{stats[4]}"
  puts "Mean coverage\t#{File.readlines(coverage)[0].to_i}"
  puts "CI\t#{benchmark_res[0]}"
  puts "accuracy\t#{benchmark_res[1]}"
  puts "ratio of detected reference\t#{(1-benchmark_res[2]).to_f.round(2)}"
end