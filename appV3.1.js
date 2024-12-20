const express = require('express');
const path = require('path');
const razers3Aligner = require('/home/tsai-apps/ProbeSearchV3.2/razers3V3'); 
const _717V2 = require('/home/tsai-apps/ProbeSearchV3.2/717V2');
const fs = require('fs');
const { exec } = require('child_process');
const { execSync } = require('child_process');
const public = __dirname + '/vue/dist/';
const app = express();
const https = require('https');


app.use(express.json());
app.use(express.static(public))

// Path to SSL certificate and private key
const privateKeyPath = '/data/private-keys/private.key';
const certificatePath = '/data/private-keys/localhost.crt';

// Read SSL certificate and private key files
const privateKey = fs.readFileSync(privateKeyPath, 'utf8');
const certificate = fs.readFileSync(certificatePath, 'utf8');
const credentials = { key: privateKey, cert: certificate };

// Create HTTPS server
const httpsServer = https.createServer(credentials, app);

app.get('/', function (req, res) {
    res.sendFile(path.join(public + "index.html"));
});

httpsServer.listen({ port: 8503 }, async () => {
    console.log('server up @ https://localhost:8503/ !')
});


/* equates user-selected db name to its path in the directory */
var db_dictionary = {
    "717V5": "/data/probesearchDB/data/717V5/g717v5_h1h2.fa",
    "PtrichocarpaV3.1": "/data/probesearchDB/data/PtrichocarpaV3.1/Ptrichocarpa_444_v3.0.fa",
    "PtrichocarpaV4.0": "/data/probesearchDB/data/PtrichocarpaV4.0/Ptrichocarpa_533_v4.0.fa",
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

const serverLogFile = path.join("/home/tsai-apps/ProbeSearchV3.2/", 'usageV3.log');

/**
 * listens for put request: calls batmis or bowtie2 child process and sends back alignment visualization
 */
app.put('/',  async (req, res) => {
	const { method, url } = req;
	
	const logEntry = `${new Date().toISOString()} - ${method} ${url}\n`;
    fs.appendFile(serverLogFile, logEntry, (err) => {
        if (err) {
            console.error('Error writing to log file:', err);
        }
    });
	//import { rAligner } from './home/webserver/probesearch/razers3.js';
    /**
     * call bowtie2 or batmap
     */
	const input = req.body.read;
	const db = req.body.db;
	const X = req.body.mismatches;
	const maxHit = req.body.maxHit;
	const gapFlag = req.body.gapFlag === "true" ? true : false;
	const aligner = req.body.aligner;
	 
	console.log("input: "+input);
	console.log("db: "+db);
	console.log("X: "+X);
	console.log("maxHit: "+maxHit);
	console.log("gapFlag: "+gapFlag);
	console.log("aligner: "+aligner);
	 
	 
    if (req.body.read.length > 50) {
        exec('./bowtie/bowtie2 -x indices/' + String(req.body.db) + '/' + String(req.body.db) + ' -k 30 -c ' + String(req.body.read) + '  --very-sensitive --no-hd' , 
            (error, stdout, stderr) => {
                if (error) {
                    console.error(`exec error: ${error}`);
                }
                // parse sam file
                //console.log(stdout);
                stdout = parse(String(req.body.read), String(req.body.db), stdout);
                res.send(stdout);
        }); 
    } else {
        let ts = Date.now();
        // make fasta input
		console.log(`Line 58: ${new Date().toLocaleTimeString()}`);
        execSync('echo ">seq\n' + String(req.body.read) + '" > ' + ts + "input.fa"); 
		
		
		if(aligner=="bowtie2"){
			
			console.log("----------In bowtie2-------------");
			if(db=="sPta717V2"){
				_717V2.DB717V2(input, X)
				.then(stdout => {
					res.send(stdout);
				})
				.catch(error => {
					console.error("Error:", error);
				});
				
			}else{
				exec('./batmis/scripts/batmap -q ./' + ts + 'input.fa -g indices/' + String(req.body.db) + '/' + String(req.body.db) + ' -n' + String(req.body.mismatches) + ' -m50 -o ' + ts + 'output.txt', 
            (error, stdout, stderr) => {
                if (error) {
                    console.error(`exec error: ${error}`);
                    execSync('rm *bin');
                }
                // remove header from sam file
                exec('grep -v ^@ ' + ts + 'output.txt', (error, stdout, stderror) => {
                    // remove input and output files 
                    execSync('rm ' + ts + 'input.fa ' + ts + 'output.txt');
                    // parse sam file
					try{
						stdout = parse(String(req.body.read), String(req.body.db), stdout);
						res.send(stdout);
					}catch(error){
						console.error('Error: ', error);
					}
					console.log(`bowtie2: ${new Date().toLocaleTimeString()}`);
                });                   
			}); 
				
			}
			
		}else{
			console.log("----------In razers3-------------");
			razers3Aligner.rAlignerV3(input, db, X, maxHit, gapFlag)
			.then(stdout => {
				res.send(stdout);
				console.log(`razers3: ${new Date().toLocaleTimeString()}`);
			})
			.catch(error => {
				console.error("Error:", error);
			});
		}
		
    }
})

/**
 * parses the SAM file input and gathers reference sequence.
 * @param sequence input read 
 * @param db database
 * @param sam SAM file  
 */
function parse(sequence, db, sam) {
    res = "" 
    sam = sam.split("\n");
	
	const prepareDataList = [];
	// Process the lines
	for (let i = 0; i < sam.length; i ++) {
		// Take the first three lines of each entry and join them with tabs
		let temp = sam[i].split('\t')
		let entry = temp[2]+'\t'+ temp[3]
		prepareDataList.push(entry);
	}
	
	let geneDict = get_gene_AtOnce(prepareDataList, sequence,db);
	let featureMapping = {"CDS": "CDS", "three_prime_UTR": "3′-UTR", "five_prime_UTR": "5′-UTR"};
	
    for (var i = 0; i < sam.length - 1; i++) {
		
        target = sam[i].split("\t");
        var strand = "+";
        var complement = false;
		let cdsFlag = false
		let cdsStart = 0
		let cdsEnd = 0
		
        if (target[2] != "*") { // valid target
            if (target[1] === "16" || target[1] === "272") { // reverse complement
                strand = "-";
                complement = true;
            }
			
            res += "\t" + target[2] + " : " + target[3] + " (" + strand + ")\n"; 
			var gene = ""
				if(!db.includes("sPta717")){
					if(Object.keys(geneDict[target[3]]).length==1){
						feature = "";
						if(geneDict[target[3]]["gene"]["complementFeatureLeft"]!=null){
							if(geneDict[target[3]]["gene"]["complementFeatureLeft"]!="")
								feature += "<br><span style='margin-left: 135px;'>"+ geneDict[target[3]]["gene"]["complementFeatureLeft"]+geneDict[target[3]]["gene"]["leftOrientation"]+"</span>";
							if(geneDict[target[3]]["gene"]["complementFeatureRight"]!="")
								feature += "<br><span style='margin-left: 135px;'>"+geneDict[target[3]]["gene"]["complementFeatureRight"]+geneDict[target[3]]["gene"]["rightOrientation"]+"</span>";
						}
						//gene = geneDict[target[3]]["gene"]["ID"] +feature;
						gene = geneDict[target[3]]["gene"]["ID"];
						
					}else{
						geneID = geneDict[target[3]]["gene"]["ID"];
						featureKeys = Object.keys(geneDict[target[3]])
						if(featureKeys.includes("CDS")){
							if(strand==="-"){
								
								cdsStart = geneDict[target[3]]["CDS"]["cdsBetweenReverseStart"]
								cdsEnd = geneDict[target[3]]["CDS"]["cdsBetweenReverseEnd"]
							}else{
								cdsStart = geneDict[target[3]]["CDS"]["cdsBetweenStart"]
								cdsEnd = geneDict[target[3]]["CDS"]["cdsBetweenEnd"]
							}
							
							feature = "<span style='background-color: chartreuse;'>"+featureMapping["CDS"]+"</span>";
							cdsFlag = true;
						}else if(featureKeys.includes("three_prime_UTR")){
							feature = featureMapping["three_prime_UTR"];
						}else if(featureKeys.includes("five_prime_UTR")){
							feature = featureMapping["five_prime_UTR"];
						}else{
							feature = "Intron";
						}
						//gene = geneID + " Feature: "+ feature;
						gene = geneID
					}
					
				}else{
					if (geneDict[target[3]] == "."){
						gene = "intergenic";
					}else{
						gene = geneDict[target[3]];
					}
				}
			
			
            res += "\tgene: " + gene + "\n";
            cigar = target[5];
            read = target[9];
            reference = String(execSync('samtools faidx ' +  db_dictionary[db] + ' ' + target[2] + ":" + target[3] 
                                        + "-" + (parseInt(target[3]) + sequence.length - 1) + ' | sed 1d')).toUpperCase();
			
            reference = reference.replace("\n", "");
            if (strand === "-") {
                read = read.split("").reverse().join("").replace("\n", ""); 
                reference = reference.split("").reverse().join("").replace("\n", "");
                read = comp(read);
                reference = comp(reference);
            }
			
            res += illustrate(cigar, read, reference, db, target[2], target[3], cdsFlag, cdsStart,cdsEnd); 
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
function illustrate(cigar, read, reference, db, chrom, pos, cdsFlag,cdsStart, cdsEnd ) {
    var ptr = 0; // alignment iterator
    var mismatches = 0; // number of mismatches 
    /* parse cigar  */
    for (var i = 0; i < cigar.length; i++) { 

        var j = cigar.length - 1;
        // check letters and grab index of closest letter
        if (cigar.substring(i, cigar.length).indexOf('M') != -1) { 
            j = cigar.substring(i, cigar.length).indexOf('M') + i;
        }
        if (cigar.substring(i, cigar.length).indexOf('I') != -1 && (cigar.substring(i, cigar.length).indexOf('I') + i) < j) { 
            j = cigar.substring(i, cigar.length).indexOf('I') + i;
        }
        if (cigar.substring(i, cigar.length).indexOf('D') != -1 && (cigar.substring(i, cigar.length).indexOf('D') + i) < j) { 
            j = cigar.substring(i, cigar.length).indexOf('D') + i;
        }
        var letter = cigar.charAt(j);
        var nucleotides = parseInt(cigar.substring(i, j));

        if (letter === 'M') {
            ptr += nucleotides;
        } else if (letter === "I") {
            var space = "-".repeat(nucleotides);
            reference = reference.substring(0, ptr) + space + reference.substring(ptr, reference.length);
            ptr += nucleotides;
        } else if (letter === "D") {
            var space = "-".repeat(nucleotides);
            read = read.substring(0, ptr) + space + read.substring(ptr, read.length);
            ptr += nucleotides;
        }
        if (j != cigar.length - 1) {
            i = j;
        } // if
    } // for
    
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
	//if(cdsFlag){
	//	// Extract the substrings before, within, and after the specified range
	//	
	//	const beforeRange = inputString.slice(0, cdsStart-1);
	//	const withinRange = inputString.slice(cdsStart-1, cdsEnd);
	//	const afterRange = inputString.slice(cdsEnd);
	//	const wrappedString = `${beforeRange}<span style="background-color: chartreuse;">${withinRange}</span>${afterRange}`;
	//	inputString = wrappedString;
	//}
	
    return illustration += "\n\t\tT   " + inputString + "\n\t" + "total mismatches: " + mismatches + "\n\n";

} // illustrate


/**
 * takes in a sequence and returns its base-pair complement 
 * @param str input read or reference 
 * @returns complement
 */
function comp(str) {
    var res = ""; 
    var complement = { "A": "T", "T": "A", "G": "C", "C": "G" }
    for (var i = 0; i < str.length; i++) {
        res += complement[str.charAt(i)]
    } // for
    return res;
} // comp


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
        var chromosome = target_info[0].split(" ")[0].replace(/\D/g, ""); 
        var mismatch = target_info[(x + 4)].substr(-2).trim();
        if (mismatch != min) {
            sorted[mismatch] = new Array(0);
            min = mismatch;
        }
        sorted[mismatch].push([chromosome, targets[i]]);
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

function get_gene_AtOnce(prepareDataList, sequence, db) {
	
	let result = {}
	const chrom = [];
	const start_coord = [];
	const end_coord = [];
	for (var i = 0; i < prepareDataList.length-1; i++) {
		target = prepareDataList[i].split("\t");
		chrom.push(target[0])
		start_coord.push(parseInt(target[1]))
		end_coord.push((parseInt(target[1]) + parseInt(sequence.length - 1)))
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
		gff3Command = 'bedtools intersect -a ' + ts + 'query.bed -b ' + featuredb_dictionary[db] + ' -wao -nonamecheck | awk \'{print $1, $2, $3, $4, $5, $6, $7, $8, $12, $13}\' > '+ ts+ 'FeatureFile.txt'
		console.log("gff3Command: "+gff3Command);
		execSync(gff3Command)
		
		complementCommand = 'bedtools intersect -a ' + ts + 'query.bed -b ' + complementFeaturedb_dictionary[db] + ' -wao -nonamecheck | awk \'{print $1, $2, $3, $4, $5, $6, $7, $8, $12, $13}\' > '+ ts+ 'ComplementFeatureFile.txt'
		console.log("complement Command: "+complementCommand);
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
				//if the gene is intergenic or not
				if(item[5]!="."){
					if(!dictFeature[key]){
						dictFeature[key] = {};
					}
					let chromID = item[8].split(";")[0].split(".")
					dictFeature[key][item[5]] = {"ID": chromID[0].replace("ID=","")+"."+chromID[1], "chromStart": item[1], "chromEnd": item[2],"featureStart": item[6], "featureEnd": item[7], "length": item[9]};
					
					if(item[5]=="CDS"){
						let cdsBetweenStart;
						let cdsBetweenEnd;
						let cdsBetweenReverseEnd;
						let cdsBetweenReverseStart;
						
						
						if(dictFeature[key][item[5]]["chromStart"]>dictFeature[key][item[5]]["featureStart"] && dictFeature[key][item[5]]["chromEnd"] < dictFeature[key][item[5]]["featureEnd"]){
							//for +ve orientation
							cdsBetweenStart = 1;
							cdsBetweenEnd = parseInt(sequence.length);
							//for -ve orientation
							cdsBetweenReverseStart = 1;
							cdsBetweenReverseEnd = parseInt(sequence.length);
						}else if(dictFeature[key][item[5]]["featureStart"]>dictFeature[key][item[5]]["chromStart"] && dictFeature[key][item[5]]["chromEnd"] < dictFeature[key][item[5]]["featureEnd"]){
							//CDS between = (featureStart-chromStart) <--> (featureEnd-chromEnd)
							//for +ve orientation
							cdsBetweenStart = parseInt(dictFeature[key][item[5]]["featureStart"] - dictFeature[key][item[5]]["chromStart"])+1;
							cdsBetweenEnd = parseInt(sequence.length); //length of chrom
							
							//for -ve orientation
							cdsBetweenReverseStart = 1;
							cdsBetweenReverseEnd = parseInt(sequence.length) - parseInt(dictFeature[key][item[5]]["featureStart"] - dictFeature[key][item[5]]["chromStart"]);
							
							
						}else if(dictFeature[key][item[5]]["chromStart"]>dictFeature[key][item[5]]["featureStart"] && dictFeature[key][item[5]]["featureEnd"]<dictFeature[key][item[5]]["chromEnd"]){
							//CDS between = (chromStart-featureStart) <--> (chromEnd-featureEnd)
							//for +ve orientation
							cdsBetweenStart = 1;
							cdsBetweenEnd = parseInt(sequence.length) - parseInt(dictFeature[key][item[5]]["chromEnd"]-dictFeature[key][item[5]]["featureEnd"]);
							
							//for -ve orientation
							cdsBetweenReverseStart = parseInt(dictFeature[key][item[5]]["chromEnd"]-dictFeature[key][item[5]]["featureEnd"])+1;
							cdsBetweenReverseEnd = parseInt(sequence.length) - parseInt(dictFeature[key][item[5]]["featureStart"] - dictFeature[key][item[5]]["chromStart"]);
							
						}else{
							//for +ve orientation
							cdsBetweenStart = parseInt(dictFeature[key][item[5]]["featureStart"] - dictFeature[key][item[5]]["chromStart"])+1;
							cdsBetweenEnd = parseInt(sequence.length) - parseInt(dictFeature[key][item[5]]["chromEnd"]-dictFeature[key][item[5]]["featureEnd"]);
							
							//for -ve orientation
							cdsBetweenReverseStart = parseInt(dictFeature[key][item[5]]["chromEnd"]-dictFeature[key][item[5]]["featureEnd"])+1;
							cdsBetweenReverseEnd = parseInt(sequence.length) - parseInt(dictFeature[key][item[5]]["featureStart"] - dictFeature[key][item[5]]["chromStart"]);
				
						}
						
						dictFeature[key][item[5]]["cdsBetweenStart"] = cdsBetweenStart;
						dictFeature[key][item[5]]["cdsBetweenEnd"] = cdsBetweenEnd;
						dictFeature[key][item[5]]["cdsBetweenReverseStart"] = cdsBetweenReverseStart;
						dictFeature[key][item[5]]["cdsBetweenReverseEnd"] = cdsBetweenReverseEnd;
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
		
	}else{
		//bedtools query .genetable file
		command = 'bedtools intersect -a ' + ts + 'query.bed -b data/' + db + '/' + db + '.gene.table -wao -nonamecheck | awk \'{print $1, $2, $3, $7, $8}\' > '+ ts+ 'Name.txt'
		execSync(command)
		// Calling the readFiles function with file name
		let arr = readFileLines(ts+ 'Name.txt');
		const dict = {}
		for (let i = 0; i < arr.length; i++) {
			key = arr[i].split(' ')[1]
			if(!dict[key]){
				dict[key] = arr[i].split(' ')[3]
			}
		}
		
		
		Object.keys(dict).forEach((key) => {
			const value = dict[key];
			});
		execSync('rm ' + ts + 'Name.txt '+ ts + 'query.bed');
		result = dict
	}
	
	
	//return [dict, dictFeature]
	return result
}