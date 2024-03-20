const { exec } = require('child_process');
const { execSync } = require('child_process');
const sqlite = require('better-sqlite3');

const public = __dirname + '/vue/717V2/dist/';

let db = new sqlite('/data/probesearchDB/data/sPta717V2.0/coverage.db');


/*********************************************************************************************************************/


/**
 * listens for put request: calls batmis or bowtie2 child process and sends back alignment visualization
 */
 
function DB717V2(input,X) {
	
	return new Promise((resolve, reject) => {
		if (input.length > 50) {
        
			alba = String(execSync('./bowtie/bowtie2 -x indices/sPta717V2.0/sPta717V2.0alba -k 30 -c ' + String(input) + '  --very-sensitive --no-hd'))
			tremula = String(execSync('./bowtie/bowtie2 -x indices/sPta717V2.0/sPta717V2.0tremula -k 30 -c ' + String(input) + '  --very-sensitive --no-hd'))
			console.log(alba)
			console.log(tremula)
			// pass parse function with both sam files 
			combined = parse(String(input), alba, tremula)
			res.send(combined)
        
		} else {

			let ts = Date.now();
			// make fasta input
			execSync('echo ">seq\n' + String(input) + '" > ' + ts + "input.fa"); 
			
			execSync('./batmis/scripts/batmap -q ./' + ts + 'input.fa -g indices/sPta717V2.0/sPta717V2.0alba -n' + String(X) + ' -m50 -o ' + ts + 'output_alba.txt')
			execSync('./batmis/scripts/batmap -q ./' + ts + 'input.fa -g indices/sPta717V2.0/sPta717V2.0tremula -n' + String(X) + ' -m50 -o ' + ts + 'output_tremula.txt')
			
			alba = String(execSync('grep -v ^@ ' + ts + 'output_alba.txt'));
			tremula = String(execSync('grep -v ^@ ' + ts + 'output_tremula.txt'));
			console.log(alba)
			
			execSync('rm ' + ts + 'input.fa ' + ts + 'output_alba.txt ' + ts + 'output_tremula.txt');
			// rm bin on error
	
			combined = parse(String(input), alba, tremula)
			//res.send(combined)
			return resolve(combined);
			
         
    }
	
	
	
	});
}
 

/**
 * parses the SAM file input and gathers reference sequence.
 * @param sequence input read 
 * @param db database
 * @param sam SAM file  
 */
function parse(sequence, sam_alba, sam_tremula) {
    res = "" 
    sam_alba = sam_alba.split("\n");
    sam_tremula = sam_tremula.split("\n");

    for (var i = 0; i < sam_alba.length - 1; i++) {
        target_alba = sam_alba[i].split("\t");
        target_tremula = sam_tremula[i].split("\t");

        var strand = "+";
        var complement = false;
        if (target_alba[2] != "*" && target_tremula[2] != "*") { // valid target
            if (target_alba[1] === "16" || target_alba[1] === "272") { // reverse complement
                strand = "-";
                complement = true;
            }
            res += "\t" + target_alba[2] + " : " + target_alba[3] + " (" + strand + ")\n"; 
            /* get gene model name */
            var gene = get_gene(target_alba[2], parseInt(target_alba[3]), (parseInt(target_alba[3]) + parseInt(sequence.length - 1)), "sPta717V2.0");
            res += "\tgene: " + gene + "\n";
            /* get CIGAR, read, and reference sequence */
            cigar_alba = target_alba[5];
            cigar_tremula = target_tremula[5];
            read_alba = target_alba[9];
            read_tremula = target_tremula[9];
            reference = String(execSync('samtools faidx /data/probesearchDB/data/sPta717V2.0/sPta717alba.fasta ' + target_alba[2] + ":" + target_alba[3] 
                                        + "-" + (parseInt(target_alba[3]) + sequence.length - 1) + ' | sed 1d')).toUpperCase();
            reference = reference.replace("\n", "");
            if (strand === "-") {
                read_alba = read_alba.split("").reverse().join("").replace("\n", ""); 
                read_tremula = read_tremula.split("").reverse().join("").replace("\n", ""); 
                reference = reference.split("").reverse().join("").replace("\n", "");
            }
            if (complement === true) {
                read_alba = comp(read_alba);
                read_tremula = comp(read_tremula);
                reference = comp(reference);
            }
            res += illustrate(cigar_alba, cigar_tremula, read_alba, read_tremula, reference, target_alba[2], target_alba[3]); 
        } else {
            return "no targets found\n\n";
        } // if
    } // for
    return sort_illustration(res);
} // parse 


/**
 * illustrates the alignment between the read & reference.
 * @param cigar CIGAR string of SAM
 * @param read sequence of read 
 * @param reference sequence of reference 
 * @returns illustration of alignment
 */
function illustrate(cigar_alba, cigar_tremula, read_alba, read_tremula, reference, chrom, pos) {
    var ptr = 0, mismatches = 0; 
    /* parse cigar  */
    for (var i = 0; i < cigar_alba.length; i++) { 
        var j = cigar_alba.length - 1;
        // check letters and grab index of closest letter
        if (cigar_alba.substring(i, cigar_alba.length).indexOf('M') != -1) { 
            j = cigar_alba.substring(i, cigar_alba.length).indexOf('M') + i;
        }
        if (cigar_alba.substring(i, cigar_alba.length).indexOf('I') != -1 && (cigar_alba.substring(i, cigar_alba.length).indexOf('I') + i) < j) { 
            j = cigar_alba.substring(i, cigar_alba.length).indexOf('I') + i;
        }
        if (cigar_alba.substring(i, cigar_alba.length).indexOf('D') != -1 && (cigar_alba.substring(i, cigar_alba.length).indexOf('D') + i) < j) { 
            j = cigar_alba.substring(i, cigar_alba.length).indexOf('D') + i;
        }
        
        var letter = cigar_alba.charAt(j), nucleotides = parseInt(cigar_alba.substring(i, j));
        if (letter === 'M') {
            ptr += nucleotides;
        } else if (letter === "I") {
            var space = "-".repeat(nucleotides);
            reference = reference.substring(0, ptr) + space + reference.substring(ptr, reference.length);
            ptr += nucleotides;
        } else if (letter === "D") {
            var space = "-".repeat(nucleotides);
            read_alba = read_alba.substring(0, ptr) + space + read_alba.substring(ptr, read_alba.length);
            ptr += nucleotides;
        }
        if (j != cigar_alba.length - 1) { i = j; } 
    } // for

    /* we redo this parsing process for tremula sam file */
    for (var k = 0; k < cigar_tremula.length; k++) { 
        var l = cigar_tremula.length - 1;
        // check letters and grab index of closest letter
        if (cigar_tremula.substring(k, cigar_tremula.length).indexOf('M') != -1) { 
            l = cigar_tremula.substring(k, cigar_tremula.length).indexOf('M') + k;
        }
        if (cigar_tremula.substring(k, cigar_tremula.length).indexOf('I') != -1 && (cigar_tremula.substring(k, cigar_tremula.length).indexOf('I') + k) < l) { 
            l = cigar_tremula.substring(k, cigar_tremula.length).indexOf('I') + k;
        }
        if (cigar_tremula.substring(k, cigar_tremula.length).indexOf('D') != -1 && (cigar_tremula.substring(k, cigar_tremula.length).indexOf('D') + k) < l) { 
            l = cigar_tremula.substring(k, cigar_tremula.length).indexOf('D') + k;
        }
        
        var letter = cigar_tremula.charAt(l), nucleotides = parseInt(cigar_tremula.substring(k, l));
        if (letter === 'M') {
            ptr += nucleotides;
        } else if (letter === "I") {
            ptr += nucleotides;
        } else if (letter === "D") {
            var space = "-".repeat(nucleotides);
            read_tremula = read_tremula.substring(0, ptr) + space + read_tremula.substring(ptr, read_tremula.length);
            ptr += nucleotides;
        }
        if (l != cigar_tremula.length - 1) { k = l; } 
    } // for

    /* build and return illustration */
    var illustration = "\t\tA   " + read_alba.substring(0, read_alba.length)
                     + "\n\t\tT   " + read_tremula.substring(0, read_tremula.length)
                     + "\n \t\t    ";
    for (var i = 0; i < read_alba.length; i++) {
        if (read_alba.substring(i, i + 1) == reference.substring(i, i + 1) && read_tremula.substring(i, i + 1) == reference.substring(i, i + 1)) { illustration += "|"; } 
        else { illustration += " "; mismatches++; } 
    } // for
    illustration += "\n\t\tQ   " + reference.substring(0, read_alba.length);
    return illustration += get_coverage(chrom, parseInt(pos), parseInt(read_tremula.length)) 
                        + "\n\t" + "total mismatches: " + mismatches + "\n\n";
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
function sort_illustration(illustration) {
    var targets = illustration.split("\n\n");
    var sorted = [];
    var min = 0;

    sorted[0] = new Array(0);
    for (var i = 0; i < targets.length - 1; i++) {
        var target_info = targets[i].split("\n");
        var chromosome = target_info[0].split(" ")[0].replace(/\D/g, ""); // grabs chromosome coordinate 
        var mismatch = target_info[7].substr(-2).trim();
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
    } // for
    return res;
} // sort_illustration


/**
 * gets the gene model name from the gff3 file via `bedtools intersect`
 * @param chrom chromosome 
 * @param start_coord starting coordinate
 * @param end_coord end coordinate  
 * @returns gene model name 
 */
 function get_gene(chrom, start_coord, end_coord, db) {
    let ts = Date.now();
    execSync('echo "' + String(chrom) + '\t' + String(start_coord) + '\t' + String(end_coord) + '" > ' + ts + "query.bed") // make input file
    var intersect = String(execSync('bedtools intersect -a ' + ts + 'query.bed -b /data/probesearchDB/data/' + db + '/' + db + '.gene.table -wb -nonamecheck | head -n 1'));
    execSync('rm ' + ts + 'query.bed'); // delete input file 

    if (intersect === "") { return "intergenic" }
    else { return intersect.split("\t")[6].replace("\n", "") }
} // get_gene


/**
 * gets the coverage score from 
 * @param chrom chromosome 
 * @param pos position
 * @param length length of read
 * @returns coverage score
 */
function get_coverage(chrom, pos, length) {
    var coverage = "\n\t\tC   ";
    let sql = `SELECT * FROM ${chrom} WHERE pos BETWEEN ${pos} AND ${(pos + (length - 1))} LIMIT ${length}`;
    if (chrom.charAt(0) === "s") { sql = `SELECT * FROM Scaffold WHERE chrom='${chrom}' AND pos BETWEEN ${pos} AND ${(pos + (length - 1))} LIMIT ${length}`; } 
    var rows = db.prepare(sql).all();

    if (rows.length != 0) {
        for (const row of rows) { 
            var score = row.score;
            if (score == 0) { coverage += "0" }
            else if (score <= 10) { coverage += "A" }
            else if (score <= 50) { coverage += "B" }
            else if (score <= 100) { coverage += "C" }
            else if (score <= 500) { coverage += "D" }
            else { coverage += "E" }
        }
    } else { coverage += "0".repeat(length); }
    return coverage;
}

module.exports = {DB717V2}; // Export the function