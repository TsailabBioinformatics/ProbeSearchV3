<script>
  export default {
    name: "SequenceForm",
    data() {
      return {
        read: '',
          db: [],
          error: ''
      }
    },
    methods: {
      /**
       * performs error checking, then emits form status to App.vue
       */
      submit() { 
        const regex = /[^ATGC]/i;
        if (this.read.search(regex) !== -1) {
          this.error = "only DNA sequences valid";
        } else if (this.read.length < 20) {
          this.error = "input a DNA sequence of length 20 or higher";
        } else if (this.db.length === 0) {
          this.error = "select a database";
        } else {
          this.error = '';
          this.$emit("valid");
        }
      }
    },
    emits: ["valid"]
  }
</script>


<template>
  
  <div class="parent">

    <div style="display: flex; justify-content: center; margin: 0;"><a href="http://aspendb.uga.edu/"><img src="../assets/aspendb_bw2.png" style="opacity: 0.9" width="100" height="60  "/></a></div>
    <div class="divider"></div> 
    <p style="width: fit-content; padding: 0 1%; border-radius: 3px; background-color: #f2f2f2">
      input primer/probe/gRNA sequence
    </p>
    <textarea rows="2" placeholder="GGGTTCTGCCAATTTAAGCCACATGGCTCAATGGGAGA" v-model="read"></textarea> 
    <p style="width: fit-content; padding: 0 1%; border-radius: 3px; background-color: #f2f2f2"> 
      select database(s)
    </p>

    <div style="display: flex; flex-direction: column;">

      <div class="checkboxes">
        <div>
          <input v-model="db" type="checkbox" name="sPta71a7V1" value="sPta717aV1">
          <label for="sPta717aV1">sPta717aV1</label>
        </div>
        <div>
          <input v-model="db" type="checkbox" name="sPta717tV1" value="sPta717tV1">
          <label for="sPta717tV1">sPta717tV1</label>
        </div>
        <div>
          <input v-model="db" type="checkbox" value="PtrichocarpaV3.1">
          <label>PtrichocarpaV3.1</label>
        </div>
        <div>
          <input v-model="db" type="checkbox" value="PtrichocarpaV4.0">
          <label>PtrichocarpaV4.0</label>
        </div>
        <div>
          <input v-model="db" type="checkbox" value="DeltoidesWV94">
          <label>DeltoidesWV94</label>
        </div>
        <div>
          <input v-model="db" type="checkbox" name="717V5" value="717V5">
          <label for="717V5">717V5</label>
        </div>
      </div>
      <button v-on:click="submit()">
        <p>search</p>
      </button>

      <!-- display on error -->
      <div class="error" v-if="error != ''">
        <div style="width: fit-content; background-color: #f5f5f7; border-radius: 2px; padding: 0 1%">
          <p> {{ error }} </p>
        </div>
      </div>

    </div>
    
  </div>
</template>


<style scoped>

p {
    margin: 1% 2% 1% 2%;
    font-size: 14px;
}
textarea {
    width: 96%;
    margin: 0 2% 0 2%;
    padding: 4%;
    font-size: 12px;
    font-weight: 500px;
    border: 1px solid var(--color-background);
    resize: none;
}
input {
    margin: 1% 1% 1% 0;
    width: fit-content;
}
button {
    width: 45%;
    margin: 4% auto 5% auto;
    color: var(--vt-c-black-mute);
    background-color: #f5f5f7; 
    border-radius: 1px;
    border: 1px solid var(--vt-c-black-mute);
    cursor: pointer;
}
button:hover {
    background-color: #f2f2f2;
}
.parent {     
    width: 90vw;                                                                                                                                             
    display: flex; 
    flex-direction: column;
    flex-shrink: 0;
    margin: 4%;
    background-color: #fbfbfb;
}
.divider {
    width: 100%; 
    border: 1px solid var(--vt-c-black-mute); 
    border-bottom-style: none; 
    margin-bottom: 2%;
}
.checkboxes {
    width: 100%;
    margin: 0 2%;
    display: flex;
    flex-direction: column;
}
.error {
    margin-top: 1%;
    display: flex;
    justify-content: center;
    align-items: center;
}

@media (min-width: 1024px) {
.parent {
    width: 70vw; 
    margin: 5% 3% 2% 3%;
    padding: 1rem 7rem 2rem 7rem;
    opacity: .95;
    box-shadow: 0px 0px 1px 0px var(--color-background); 
    border-radius: 3px;
}
.checkboxes {
    margin: 0;
}
p {
    margin: 1% 0 1% 0;
    font-size: 18px;
}
textarea {
    width: 100%;
    padding: 4px;
    margin: 0;
    background-color: rgba(247, 247, 247, 0.90);
    font-size: 16px;
    font-weight: 500px;
}
input {
    margin: 0.5%;
    width: fit-content;
}
label {
    color: var(--vt-c-black-mute);
    font-size: 16px;
}
button {
    width: 25%;
    padding: 5px;
    margin: 1% 0;
    font-size: 14px;
    margin-left: auto;
    margin-right: auto;
} }

@media (min-width: 2048px) {
.divider {
    margin-bottom: 1%;
}
textarea {
    width: 1000px;
}
button {
    width: 15%;
}
input {
    margin: 0.5%;
} }

</style>