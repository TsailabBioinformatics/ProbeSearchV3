const fs = require('fs');
const { exec } = require('child_process');
const { execSync } = require('child_process');


/* equates user-selected db name to its path in the directory */
var db_dictionary = {
    "717V5": "/data/probesearchDB/data/717V5/g717v5_h1h2.fa",
    "PtrichocarpaV3.1": "/data/probesearchDB/data/PtrichocarpaV3.1/Ptrichocarpa_444_v3.0.fa",
    "PtrichocarpaV4.0": "/data/probesearchDB/data/Ptrichocarpa_533_v4.0.fa",
    "DeltoidesWV94": "/data/probesearchDB/data/DeltoidesWV94/PdeltoidesWV94_445_v2.0.fa",
    "sPta717V1" : "/data/probesearchDB/data/sPta717V1/Pta717s_v1.1.fa",
    "sPta717tV1" : "/data/probesearchDB/data/sPta717tV1/sPta717tremula.fasta",
    "sPta717aV1" : "/data/probesearchDB/data/sPta717aV1/sPta717alba.fasta",
    "sPta717V2.0" : "/data/probesearchDB/data/sPta717V2.0/sPta717_v2.0.fa",
    "717V5M147masked" : "/data/probesearchDB/data/717V5M147masked/g717v5_h1h2.M147round.hardmasking.fa"
};

var featuredb_dictionary = {
    "717V5": "/data/probesearchDB/data/717V5/gv5.h1.gene_exon.renamed.gff3",
    "PtrichocarpaV3.1": "/data/probesearchDB/data/PtrichocarpaV3.1/Ptrichocarpa_444_v3.1.gene_exons.gff3",
    "PtrichocarpaV4.0": "/data/probesearchDB/data/PtrichocarpaV4.0/Ptrichocarpa_533_v4.1.gene_exons.gff3",
    "DeltoidesWV94": "/data/probesearchDB/data/DeltoidesWV94/PdeltoidesWV94_445_v2.1.gene_exons.gff3",
    "717V5M147masked" : "/data/probesearchDB/data/717V5/gv5.h1.gene_exon.renamed.gff3"
};

var complementFeaturedb_dictionary = {
    "717V5": "/data/probesearchDB/data/717V5/g717v5_h1h2_ComplementDB.bed",
    "PtrichocarpaV3.1": "/data/probesearchDB/data/PtrichocarpaV3.1/Ptrichocarpa_444_v3.1_ComplementDB.bed",
    "PtrichocarpaV4.0": "/data/probesearchDB/data/PtrichocarpaV4.0/Ptrichocarpa_533_v4.1_ComplementDB.bed",
    "DeltoidesWV94": "/data/probesearchDB/data/DeltoidesWV94/PdeltoidesWV94_445_v2.1_ComplementDB.bed",
    "717V5M147masked" : "/data/probesearchDB/data/717V5/g717v5_h1h2_ComplementDB.bed"
};

/*********************************************************************************************************************/


/**
 * listens for put request: calls batmis or bowtie2 child process and sends back alignment visualization
 */
 
//app.put('/',  async (req, res) => {
//app.put('/',  async (req, res) => {
function rAlignerV3(input,db,X,maxHit,gapFlag) {
	
	return new Promise((resolve, reject) => { 
		
		
		if (input.length > 50) {
			exec('./bowtie/bowtie2 -x indices/' + String(req.body.db) + '/' + String(req.body.db) + ' -k 30 -c ' + String(req.body.read) + '  --very-sensitive --no-hd' , 
				(error, stdout, stderr) => {
					if (error) {
						console.error(`exec error: ${error}`);
					}
					// parse sam file
					console.log(stdout);
					stdout = parse(String(req.body.read), String(req.body.db), stdout);
					res.send(stdout);
			}); 
		} 
		else {
				let ts = Date.now();
				let gap = "";
				if(!gapFlag){
					gap = "-ng";
				}
				
				let errorPercentage = Math.floor((1-(X/input.length) )*100 )
				//console.log(`Start --> Line 56: ${new Date().toLocaleTimeString()}`);
				// make fasta input
				execSync('echo ">seq\n' + input + '" > ' + ts + "input.fa"); 
				exec('razers3 -i '+ errorPercentage +' -m '+ maxHit +' -tc 4 --a '+ gap +' -o '+ ts + 'output.razers /data/probesearchDB/data/'+String(db)+'/'+'*.fa ' + ts +'input.fa', 
					(error, stdout, stderr) => {
						//console.log(`After syntax complete --> Line 87: ${new Date().toLocaleTimeString()}`);
						if (error) {
							console.error(`exec error: ${error}`);
							execSync('rm *bin');
						}
						// remove header from razers file
						exec('grep -v ^@ ' + ts + 'output.razers', (error, stdout, stderror) => {
							// remove input and output files 
							execSync('rm ' + ts + 'input.fa ' + ts + 'output.razers');
							// parse razers file
							stdout = parse(String(input), String(db), stdout);
							//res.send(stdout);
							//console.log(`Line 77: ${new Date().toLocaleTimeString()}`);
							return resolve(stdout);
							//resolve(result);
						});                   
				}); 
				
			}
	
	});

}


/**
 * parses the razers file input and gathers reference sequence.
 * @param sequence input read 
 * @param db database
 * @param razers razers file  
 */
function parse(sequence, db, razers) {
    res = "" 
    razers = razers.split("\n");
	
	const prepareDataList = [];
	// Process the lines
	for (let i = 0; i < razers.length; i += 3) {
		// Take the first three lines of each entry and join them with tabs
		const entry = razers.slice(i, i + 3).join('\t');
		prepareDataList.push(entry);
	}

	geneDict = get_gene_AtOnce(prepareDataList, sequence,db);
	let featureMapping = {"CDS": "CDS", "three_prime_UTR": "3′-UTR", "five_prime_UTR": "5′-UTR"};
	
    for (var i = 0; i < prepareDataList.length - 1; i++) {
        target = prepareDataList[i].split("\t");
		
        var strand = "+";
		
		var complement = false;
		let cdsFlag = false
		let cdsStart = 0
		let cdsEnd = 0
		
		//Check if the target is Reverse or Forward
        if (target.length != 0) { // valid target
		
			if(target[3] === "R"){
				strand = "-";
			}
            
            res += "\t" + target[4] + " : " + (parseInt(target[5])+1) + " (" + strand + ")\n";
			var gene = "";
			
			//try{
				if(!db.includes("sPta717")){
					
					//---------------------------------------------------------
					if(Object.keys(geneDict[(parseInt(target[5])+1)]).length==1){
						feature = "";
						if(geneDict[(parseInt(target[5])+1)]["gene"]["complementFeatureLeft"]!=null){
							if(geneDict[(parseInt(target[5])+1)]["gene"]["complementFeatureLeft"]!="")
								feature += "<br><span style='margin-left: 135px;'>"+ geneDict[(parseInt(target[5])+1)]["gene"]["complementFeatureLeft" ]+geneDict[(parseInt(target[5])+1)]["gene"]["leftOrientation"]+"</span>"
							if(geneDict[(parseInt(target[5])+1)]["gene"]["complementFeatureRight"]!="")
								feature = "<br><span style='margin-left: 135px;'>"+geneDict[(parseInt(target[5])+1)]["gene"]["complementFeatureRight"]+geneDict[(parseInt(target[5])+1)]["gene"]["rightOrientation"]+"</span>";
						}
						//gene = geneDict[(parseInt(target[5])+1)]["gene"]["ID"] + feature;
						gene = geneDict[(parseInt(target[5])+1)]["gene"]["ID"];
						
						
						
					}else{
						geneID = geneDict[(parseInt(target[5])+1)]["gene"]["ID"];
						featureKeys = Object.keys(geneDict[(parseInt(target[5])+1)])

						if(featureKeys.includes("CDS")){
							cdsStart = geneDict[(parseInt(target[5])+1)]["CDS"]["cdsBetweenStart"]
							cdsEnd = geneDict[(parseInt(target[5])+1)]["CDS"]["cdsBetweenEnd"]
							feature = "<span style='background-color: chartreuse;'>"+featureMapping["CDS"]+"</span>";
							cdsFlag = true;
						}else if(featureKeys.includes("three_prime_UTR")){
							feature = featureMapping["three_prime_UTR"];
						}else if(featureKeys.includes("five_prime_UTR")){
							feature = featureMapping["five_prime_UTR"];
						}else{
							feature = "Intron";
						}
						gene = geneID + " Feature: "+ feature;
						//gene = geneID;
					}
					
					//---------------------------------------------------------
				}else{
					if (geneDict[(parseInt(target[5])+1)] == "."){
						gene = "intergenic";
					} else{
						gene = geneDict[(parseInt(target[5])+1)];
					}
					
				}
				
			//}catch(error){
			//	continue;
			//}
			
            res += "\tgene: " + gene + "\n";
			read = target[8].split("#Read:   ")[1]
			reference = target[9].split("#Genome: ")[1]
			
			
			res += illustrate(read, reference, cdsFlag, cdsStart,cdsEnd); 
			
			
			
        } else {
            res += "no targets found\n\n";
            break;
        } // if
    } // for
	
    res = sort_illustration(res, 1); 
    return res;
} // parse 

/**
 * illustrates the alignment between the read & reference.
 * @param cigar CIGAR string of SAM
 * @param read sequence of read 
 * @param reference sequence of reference 
 * @returns illustration of alignment
 */
function illustrate(read, reference, cdsFlag,cdsStart, cdsEnd) {

    var mismatches = 0; // number of mismatches 
    /* build illustration  */
    var min = Math.min(read.length, reference.length);
    var illustration = "\t\tQ   " + read.substring(0, min) + "\n \t\t    "; 
    for (var i = 0; i < min; i++) {
        if (read.substring(i, i + 1) == reference.substring(i, i + 1)) {
            illustration += "|";
        } else {
            illustration += " ";
            mismatches++;
        } // if
    } // for
	
	
	
	inputString = reference.substring(0, min);
	if(cdsFlag){
		// Extract the substrings before, within, and after the specified range
		
		const beforeRange = inputString.slice(0, cdsStart-1);
		const withinRange = inputString.slice(cdsStart-1, cdsEnd);
		const afterRange = inputString.slice(cdsEnd);
		const wrappedString = `${beforeRange}<span style="background-color: chartreuse;">${withinRange}</span>${afterRange}`;
		//const wrappedString = `${beforeRange}<strong style="background-color: chartreuse;">${withinRange}</strong>${afterRange}`;
		inputString = wrappedString;
	}
	
	
    //return illustration += "\n\t\tT   " + reference.substring(0, min) + "\n\t" + "total mismatches: " + mismatches + "\n\n";
    return illustration += "\n\t\tT   " + inputString + "\n\t" + "total mismatches: " + mismatches + "\n\n";

} // illustrate


/**
 * takes in the illustration created from standard output,
 * and sorts its components chronologically by mismatch, then chromosome.
 * @param illustration 
 * @returns sorted illustration
 */
function sort_illustration(illustration, x) {
    var targets = illustration.split("\n\n");
    var sorted = [];
    var min = 0;
	
    sorted[0] = new Array(0);
    for (var i = 0; i < targets.length - 1; i++) {
        var target_info = targets[i].split("\n");
		if(target_info.length == 6){
			var chromosome = target_info[0].split(" ")[0].replace(/\D/g, ""); 
			var mismatch = target_info[5].substr(-2).trim();
			
			try {
				sorted[mismatch].push([chromosome, targets[i]]);
			} catch (error) {
				sorted[mismatch] = new Array(0);
				sorted[mismatch].push([chromosome, targets[i]]);
			}
		}
    } // for
    var res = "";
    var count = 0;
    for (var j = 0; j < sorted.length; j++) {
        if (sorted[j]) { 
            sorted[j] = sorted[j].sort(); 
            for (var k = 0; k < sorted[j].length; k++) {
                res += "hit " + (++count) + ")\n" + sorted[j][k][1] + "\n\n"
            }
        }
    }
    return res;
}


/**
 * gets the gene model name from the gff3 file via `bedtools intersect`
 * @param chrom chromosome 
 * @param coord starting coordinate 
 * @returns gene model name 
*/
function get_gene(chrom, start_coord, end_coord, db) {
    let ts = Date.now();
    execSync('echo "' + String(chrom) + '\t' + String(start_coord) + '\t' + String(end_coord) + '" > ' + ts + "query.bed") 

    var intersect = String(execSync('bedtools intersect -a ' + ts + 'query.bed -b data/' + db + '/' + db + '.gene.table -wb -nonamecheck | head -n 1'));
    execSync('rm ' + ts + 'query.bed');

    if (intersect === "") { return "intergenic" }
    else { 
        return intersect.split("\t")[6].replace("\n", "") 
    }
} // get_gene 



function get_gene_AtOnce(prepareDataList, sequence, db) {
	
	let result = {}
	const chrom = [];
	const start_coord = [];
	const end_coord = [];
	for (var i = 0; i < prepareDataList.length-1; i++) {
		target = prepareDataList[i].split("\t");
		
		chrom.push(target[4])
		start_coord.push(parseInt(target[5])+1)
		end_coord.push((parseInt(target[5])+1 + parseInt(sequence.length - 1)))
	}
	
	let ts = Date.now();
	let listOfChrom = ""
	for (var i = 0; i < chrom.length; i++) {
		listOfChrom += String(chrom[i]) + '\t' + String(start_coord[i]) + '\t' + String(end_coord[i])+'\n';
		
	}
	
    execSync('cat <<EOL > '+ts + 'query.bed \n'+ listOfChrom + ' \nEOL');
	
	const readFileLines = filename =>
		fs.readFileSync(filename)
		.toString('UTF8')
		.split('\n');
		
	//bedtools query .gff3 file if(db!="sPta717")
	if(!db.includes("sPta717")){
		//gff3Command = 'bedtools intersect -a ' + ts + 'query.bed -b data/' + db + '/Ptrichocarpa_533_v4.1.gene_exons.gff3 -wao -nonamecheck | awk \'{print $1, $2, $3, $4, $5, $6, $7, $8, $12}\' > '+ ts+ 'FeatureFile.txt'
		gff3Command = 'bedtools intersect -a ' + ts + 'query.bed -b ' + featuredb_dictionary[db] + ' -wao -nonamecheck | awk \'{print $1, $2, $3, $4, $5, $6, $7, $8, $12, $13}\' > '+ ts+ 'FeatureFile.txt'
		console.log("----gff3Command----"+gff3Command);
		execSync(gff3Command)
		
		complementCommand = 'bedtools intersect -a ' + ts + 'query.bed -b ' + complementFeaturedb_dictionary[db] + ' -wao -nonamecheck | awk \'{print $1, $2, $3, $4, $5, $6, $7, $8, $12, $13}\' > '+ ts+ 'ComplementFeatureFile.txt'
		console.log("----complementCommand----"+complementCommand);
		execSync(complementCommand)
		
		//Complement Feature 
		const dictComplementFeature = {}	
		let arrComplementFeature = readFileLines(ts+ 'ComplementFeatureFile.txt');		
		for (let i = 0; i < arrComplementFeature.length; i++) {
			item = arrComplementFeature[i].split(' ');
			if(item[6]!="."){
				dictComplementFeature[item[1]] = {};
				dictComplementFeature[item[1]]["start"] = item[1]
				dictComplementFeature[item[1]]["end"] = item[2]
				dictComplementFeature[item[1]]["pos1"] = item[4]
				dictComplementFeature[item[1]]["pos2"] = item[5]
				dictComplementFeature[item[1]]["id"] = item[6]
				dictComplementFeature[item[1]]["length"] = item[7]
				
			}
			
		}
		
		//Regular gene Feature 
		const dictFeature = {}	
		let arrFeature = readFileLines(ts+ 'FeatureFile.txt');		
		for (let i = 0; i < arrFeature.length; i++) {
			item = arrFeature[i].split(' ');
			if(item.length >1){
				let key = item[1];
				if(item[5]!="."){
					if(!dictFeature[key]){
						dictFeature[key] = {};
					}
					let chromID = item[8].split(";")[0].split(".")
					dictFeature[key][item[5]] = {"ID": chromID[0].replace("ID=","")+"."+chromID[1], "chromStart": item[1], "chromEnd": item[2],"featureStart": item[6], "featureEnd": item[7], "length": item[9]};
					
					if(item[5]=="CDS"){
						let cdsBetweenStart;
						let cdsBetweenEnd;
						if(dictFeature[key][item[5]]["chromStart"]>dictFeature[key][item[5]]["featureStart"] && dictFeature[key][item[5]]["chromEnd"] < dictFeature[key][item[5]]["featureEnd"]){
							cdsBetweenStart = 1;
							cdsBetweenEnd = parseInt(sequence.length);
						}else if(dictFeature[key][item[5]]["featureStart"]>dictFeature[key][item[5]]["chromStart"] && dictFeature[key][item[5]]["chromEnd"] < dictFeature[key][item[5]]["featureEnd"]){
							//CDS between = (featureStart-chromStart) <--> (featureEnd-chromEnd)
							cdsBetweenStart = (dictFeature[key][item[5]]["featureStart"] - dictFeature[key][item[5]]["chromStart"])+1;
							cdsBetweenEnd = parseInt(sequence.length); //length of chrom
						}else if(dictFeature[key][item[5]]["chromStart"]>dictFeature[key][item[5]]["featureStart"] && dictFeature[key][item[5]]["featureEnd"]<dictFeature[key][item[5]]["chromEnd"]){
							//CDS between = (chromStart-featureStart) <--> (chromEnd-featureEnd)
							cdsBetweenStart = 1;
							cdsBetweenEnd = (dictFeature[key][item[5]]["chromEnd"]-dictFeature[key][item[5]]["featureEnd"])+1;
						}else{
							cdsBetweenStart = (dictFeature[key][item[5]]["featureStart"] - dictFeature[key][item[5]]["chromStart"])+1;
							cdsBetweenEnd = (dictFeature[key][item[5]]["chromEnd"]-dictFeature[key][item[5]]["featureEnd"])+1;
				
						}
						dictFeature[key][item[5]]["cdsBetweenStart"] = cdsBetweenStart;
						dictFeature[key][item[5]]["cdsBetweenEnd"] = cdsBetweenEnd;
					}
					
				}else{
					dictFeature[key] = {}
					if(dictComplementFeature[key]!= null){
						leftdiff = dictComplementFeature[key]["start"] - dictComplementFeature[key]["pos1"];
						rightdiff = dictComplementFeature[key]["pos2"] - dictComplementFeature[key]["end"];
						leftInfo = "";
						rightInfo = "";
						leftOrientation = "";
						rightOrientation = "";
						temp = dictComplementFeature[key]["id"].split(";")
						
						if(temp[0].split("|")[0]!="."){
							if(leftdiff>=1000){
								leftInfo = ">1000 bp to Upstream gene: " + temp[0].split("|")[0]
								leftOrientation = "("+temp[0].split("|")[1]+")";
							}else{
								if(temp[0].split("|")[1]=="-"){
									leftInfo = leftdiff+" bp to Upstream gene: <span style='background-color: Aqua;'>" + temp[0].split("|")[0]+"</span>";
									leftOrientation = "("+temp[0].split("|")[1]+")<span style='color: orangered; font-size: 12px;'>*Possible promoter</span>";
								}else{
									leftInfo = leftdiff+" bp to Upstream gene: " + temp[0].split("|")[0]
									leftOrientation = "("+temp[0].split("|")[1]+")";
								}
							}	
							
						}
						
						if(temp[1].split("|")[0]!="."){
							if(rightdiff>=1000){
								rightInfo = ">1000 bp to Downstream gene: " + temp[1].split("|")[0]
								rightOrientation = "("+temp[1].split("|")[1]+")";
							}else{
								if(temp[1].split("|")[1]=="+"){
									rightInfo = rightdiff+" bp to Downstream gene: <span style='background-color: yellow;'>" + temp[1].split("|")[0]+"</span>";
									rightOrientation = "("+temp[1].split("|")[1]+")<span style='color: orangered; font-size: 12px;'>*Possible promoter</span>";
								}else{
									rightInfo = rightdiff+" bp to Downstream gene: " + temp[1].split("|")[0]
									rightOrientation = "("+temp[1].split("|")[1]+")";
								}
							}
						}
						
						if(temp[0].split("|")[0]=="." && temp[1].split("|")[0]=="."){
							dictFeature[key]["gene"] = {"ID": "intergenic"}
						}else{
							dictFeature[key]["gene"] = {"ID": "intergenic", "complementFeatureLeft": leftInfo, "complementFeatureRight": rightInfo, "leftOrientation": leftOrientation, "rightOrientation": rightOrientation};
						}
						
					}else{
						dictFeature[key]["gene"] = {"ID": "intergenic"}
					}
					
				}
				
			}
		}
			
		execSync('rm ' + ts + 'FeatureFile.txt '+ ts + 'query.bed ' + ts + 'ComplementFeatureFile.txt');
		result = dictFeature
	}
	
	
	
	return result
}


module.exports = {rAlignerV3}; // Export the function



