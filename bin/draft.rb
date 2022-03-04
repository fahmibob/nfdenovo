def divide_pair_into_n(globListData, nAmount, fwdIdentifier)
  allData=Dir.glob(globListData)
  pair_ch=[]
  allData.each do |x|
    if x.match("_#{fwdIdentifier}")
      pair_ch.push(x.match("_#{fwdIdentifier}").pre_match.split("/")[-1])
    end
  end
  pair_ch=pair_ch.shuffle
  equalarray=pair_ch.each_slice(nAmount).to_a.count

  chunkArray=pair_ch.each_slice(equalarray).to_a
  for i in 0..chunkArray.count-1
    write=File.open("#{File.dirname(__FILE__)}/list#{i+1}.txt", "w")
    for k in 0..chunkArray[i].count-1
      write.puts(chunkArray[i][k])
    end
  end
end
