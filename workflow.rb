require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/matrix'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/ROMA'

module ROMA
  extend Workflow

  def self.run_java(data_file, module_file, output_dir, options = {})

    bar = options.delete :bar

    options = {} if options.nil?
    options["-dataFile"] = data_file
    options["-moduleFile"] = module_file
    options["-outputFolder"] = output_dir

    cp_str = [Rbbt.jar["VDAOEngine.jar"].find, Rbbt.jar["ROMA.jar"].find] * ":"
    cmd = "java -cp #{cp_str} fr.curie.ROMA.ModuleActivityAnalysis"

    CMD.cmd(cmd, options.merge(:log => true))
  end

  def self.run_R(data_file, module_file, output_dir, options = {})
    require 'rbbt/util/R'

    bar = options.delete :bar

    modules = {}
    TSV.traverse module_file, :type => :array do |line|
      next if line =~ /^#/
      name, desc, *gene_strs = line.split("\t")

      genes = []
      weights = []
      gene_strs.each do |gene|
        if m = gene.match(/(.*)\[(.*)\]/)
          gene, weight = m.values_at 1, 2
        else
          weight = "NA"
        end
        genes << gene
        weights << weight
      end
      modules[name] = {"Name" => name, "Description" => desc, "Genes" => genes, "Weights" => weights}
    end

    m = {}
    modules.keys.each do |k| m[k] = modules[k] end
    modules = m

    TmpFile.with_file(modules.to_json, false) do |module_json|
      script = <<-EOF
library(rRoma)
library(jsonlite)

data = rbbt.tsv('#{data_file}', comment.char="")
modules = fromJSON('#{module_json}')

Data.FC <- rRoma.R(ExpressionMatrix = data, ModuleList = modules, FixedCenter = TRUE, ExpFilter=FALSE, PCSignMode="CorrelateAllWeightsByGene", GeneOutDetection='L1OutSdMean')

selected = SelectGeneSets(Data.FC)

matrix = Data.FC$SampleMatrix[selected,]

data = as.data.frame(t(matrix))

       EOF

      if bar
        init = false
        monitor = Proc.new do |line| 
          if line =~ /\[\d+\/(\d+)\] Working on/ 
            bar.max = $1.to_i unless init
            init = true
            bar.tick 
          end
        end
      else
        monitor = false
      end
      matrix = R.run script, :monitor => monitor

      matrix
    end
  end

  def self.run(*args)
    run_java(*args)
  end

  helper :file_path do |content, input|
    filename = file('inputs')[input].find
    Open.write(filename, content.to_s)
    filename
  end

  input :data, :text, "Data"
  task :data_matrix => :tsv do |data|
    tsv = TSV === data ? data : TSV.open(data)
    tsv = tsv.to_list{|v| Misc.mean(v.collect{|v| v.to_f}) } if tsv.type == :double

    stream = StringIO.new
    stream << tsv.all_fields * "\t" << "\n"
    TSV.traverse tsv, :into => stream do |k,vs|
      next if k.nil? or k.empty?
      k + "\t" << vs * "\t" << "\n"
    end
    stream.rewind
    stream
  end

  dep :data_matrix
  input :modules, :tsv, "Modules"
  task :run => :array do |modules|
    data_file = step(:data_matrix).path
    module_file = file_path modules, :modules

    output_dir = file('output')
    bar = self.progress_bar :desc => "Performing ROMA"
    ROMA.run(data_file, module_file, output_dir, :bar => bar)
    
    output_dir.glob "*"
  end

  dep :run
  task :profile => :tsv do
    lines = step(:run).file('output/moduletable_simple_T.dat').read.split("\n")
    headers, content = lines.shift.split(" ")
    lines.shift
    fields = []
    (headers.to_i-1).times do 
      fields << lines.shift.split("\t").first
    end
    tsv = TSV.setup({}, :key_field => "Sample", :fields => fields, :type => :list, :cast => :to_f)
    while line = lines.shift
      key, *values = line.split("\t")
      tsv[key] = values[0..fields.length-1]
    end

    tsv
  end

  input :data, :text, "Data"
  task :normalize => :tsv do |data|
    tsv = TSV.open(data).change_key("Associated Gene Name", :identifiers => Organism.identifiers(Organism.default_code("Hsa"))).to_list{|v| Misc.mean(v.collect{|v| v.to_f}) }
    fields = tsv.fields

    if fields.length != fields.uniq.length

      new = tsv.annotate({})
      new.type = :double
      new.fields = fields.uniq

      tsv.each do |k,vs|
        nvsh = {}
        vs.zip(fields).each do |v,f|
          nvsh[f] ||= []
          nvsh[f] << v
        end

        nvs = new.fields.collect do |f| 
          Misc.mean(nvsh[f])
        end

        new[k] = nvs
      end

      tsv = new
    end
    tsv.to_s
  end
end

#require 'ROMA/tasks/basic.rb'

#require 'rbbt/knowledge_base/ROMA'
#require 'rbbt/entity/ROMA'

