# ProbeSearchV3
[//]: <TODO: figure out the best way to structure and tranfer `indices/` and `data/` directories>

#### Overview
Probesaerch is a tool designed to facilitate the alignment of primer/probe/gRNA sequences against selected genome databases. It currently has two versions: V3 and V3.2. The app is written in Nodejs, Expressjs, and Vue.js. 

#### Versions
1. **V3**: 
   - Allows users to input primer/probe/gRNA sequences and select a set of genome databases.
   - Upon submission, the app returns the result of the input read aligned against the selected genome using the Bowtie2 software.
   - Provides an option to change the aligner to Razers3, primarily used for aligning with short reads.

2. **V3.2**:
   - Includes all features of V3.
   - Additionally provides Gene Features and CDS Overlapping Info.

### Server Information
The app is currently hosted on server 172.30.18.104 with the following configurations:
- V3: Running on port 8503, accessible through http://tsailab.gene.uga.edu/probesearch/v3/
- V3.2: Running on port 8090, accessible through http://tsailab.gene.uga.edu/probesearch/v3.2/

### Requirements
- Node.js: https://nodejs.dev/learn/how-to-install-nodejs

#### Installation and Execution of Razers3
- Follow the instructions provided in the [SeqAn GitHub repository](https://github.com/seqan/seqan/tree/main/apps/razers3) to set up Razers3 aligner executable.
  
### Instructions
1. Clone this repository. \
`git clone https://github.com/TsailabBioinformatics/ProbeSearchV3`

2. Travel to the parent project directory and install the required node packages. \
`npm install` \
Then change into the frontend `vue/` directory, and install its packages. \
`cd vue/` \
`npm install`

3. Make the directories `data/` and `indices/` directories to the parent directory. Refer to [Editing Backend section.](#editting-probesearchv3)

4. Run the app. \
`node app.js` 

---

### Understanding ProbeSearchV3
Here is description of full stack of ProbeSearchV3. 

**Directory Structure**
- `data/`: Fasta files are located in `/data/probesearchDB/data/` on the server
- `bowtie/`: Bowtie2 executables are located in `bowtie/`  
- `scripts/`: Helper scripts for Bowtie2
- `indices/`: Bowtie2 uses indexed versions of the fasta files. Indices are located in `indices/`. 
  - In order to create an index for bowtie, use `./bowtie/bowtie2-build <path_to_fasta> <index_name>`. Then move these indices to `indices/`
  - In order to create an index for batmis, use `./usr/local/bin/build_index <path_to_fasta> <index_name>`. Then move these indices to `indices/`
  - Indices are integral to bowtie2 alignment. `app.js` makes a call to `./bowtie/bowtie2 -x indices/<db> -k 30 -c <read>`. Here, `<db>` represents the user-selected database, and `<read>` represents the input read.
  - `GeneTable/`: contains gene database used by V3.2 
- `vue/`: Frontend Vue files

**Data Flow** \
The following details how data moves in the app. We'll see how the client-side form submission *request* leads to the eventual alignment result *response*. 

1. Initial Input Processing \
  User input is first handled by `SequenceForm.vue`. Within `SequenceForm.vue`, these data are tracked as `read: ''` and `db: []`. `read` is a string because it represents a primer/probe/gRNA sequence, and `db` is an array because it represents a set of genome databases. When the user fills out the form, the variables are defined respetively. Then, when the user submits the form (e.g., clicks `Search`), the variables are propagated, or *emitted* (in Vue lingo), up to its parent component `App.vue`. `App.vue`, then processes the user input and makes a put request to the Express api, `app.js` (`axios.put('/', payload)`). 
  
2. API Functions \
  The put function first calls a child process: `./bowtie/bowtie2 -x indices/ + String(req.body.db) -k 30 -c + String(req.body.read) + --end-to-end --no-hd`. Note, `req.body.db` and `req.body.read` is the Express way of accessing the inputted db and read respectively. This child process initiates the bowtie2 alignment function, which returns a SAM file as `stdout`. Then, the SAM file, or `stdout` is parsed and molded into an alignment illustration. The illustration of the alignmened is thne returned to `App.vue`. 
  
3. Returning Alignment Result \
  `app.js` sends back data in the form of an illustration to `App.vue`. In turn, `App.vue` passes this data to a child component called `AlignmentResult.vue`. Once the `AlignmentResult` component(s) receive the data, they show on client side. 

#### UI
Both versions share the same user interface to avoid complexity. However, their backend implementations differ.

#### Scripts
Each version has three separate scripts:
1. **Bowtie2 Scripts**:
   - Handles alignment using Bowtie2 for all genome databases except sPta717V2.0.
   - V3: `app.js`
   - V3.2: `appV3.js`

2. **Specialized Bowtie2 Script**:
   - Specifically for the sPta717V2.0 database due to differences in output format.
   - Script: `717V2.js`

3. **Razers3 Scripts**:
   - Utilizes Razers3 aligner for alignment.
   - V3: `razers3.js`
   - V3.2: `razers3V3.js`

#### Entry Points
- V3: Entry point is through `app.js`.
- V3.2: Entry point is through `appV3.js`.
- 

#### Editting ProbeSearchV3
- **Frontend**
  - The frontend can be editted from within the `vue/` directory. 
  - The three components as of now are the parent, `App.vue`, and its children, `SequenceForm.vue` & `AlignmentResult.vue`. 
  - After editting the vue files, run `npm run build` from within the `vue/` directory. This command builds the static HTML, css, and js for the frontend and places them in `vue/dist/` directory. The express app, `app.js`, references these static files automatically, so that's all: edit then build. 
- **API**
  - The main app can be editted at `app.js` or `appV3.js` based on the version of the app. It is written in Express.
- **Backend**
  - The backend consists of two mains parts: the fasta files & the indexed fasta files. In order to edit the backend, or implement more genomes, do the following:
  1. Transfer whatever fasta files you wish to the parent directory.
  2. Build the indices: `./bowtie/bowtie2-build <path_to_fasta> <index_name>`. Name the index (second parameter) the same name as you display the database on the frontend.
  3. Move the raw fasta file to the `data/` directory, and move the indexed fasta files to the `indices/` directory.
  4. Go the Express `app.js`, and add to `db_dictionary`, which facilitates a dictionary between the frontend name of the database and the location of the fasta file.
