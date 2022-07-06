<script>
  export default {
    name: "AlignmentResult",
    props: ['header', 'data', 'id', 'total'],
    methods: {
      swapleft() {this.$emit("swap-left") },
      swapright() { this.$emit("swap-right") },
      reload() { 
        window.location.reload();
        window.pageYOffset || document.documentElement.scrollTop || document.body.scrollTop;
      },
      expand() { this.$emit("expand") },
      download() { /* TODO */ }
    },
    emits: ["swap-left", "swap-right", "expand"]
  }
</script>


<template>
  <div class="parent">

    <p class="expand" @click="expand()"> ⛶ </p>
    <p class="instruction"> maximize </p>
    <div class="header"> 
        <div>
          <p class="swap-left" v-if="id > 1" @click="swapleft()"> &lt; </p> <p v-else> &nbsp; </p> 
          <p class="instruction"> previous result </p>
        </div>
        <h3 v-if="this.total > 1">Alignment Result {{id}} </h3> 
        <div>
          <p class="swap-right" v-if="id < total" @click="swapright()"> > </p> <p v-else> &nbsp;</p>
          <p class="instruction"> next result </p>
        </div>
    </div>

    <pre style="white-space: pre-line;">
        <h3><b>{{header}}</b></h3>
    </pre>
    <div class="divider"></div> 
    <pre>
        <p>{{data}}</p>
    </pre>

    <div class="footer">
        <!-- TODO <p class="download" @click="download"> download </p> -->
        <p class="search-again" @click="reload"> search again </p>
    </div>
    
  </div>
</template>


<style scoped>

.parent {
    width: 90vw;  
    display: flex; 
    flex-direction: column;
    flex-shrink: 0;
    margin: 0;
    background-color: #fbfbfb;
    border: 1px solid var(--vt-c-black-mute);
}
.header {
    display: flex;
    justify-content: space-evenly;
    align-items: center;
}
.divider {
    width: 100%; 
    border: 1px solid var(--vt-c-black-mute); 
    border-bottom-style: none; 
    margin: 0;
}
.footer {
    display: flex;
    align-items: flex-end;
    justify-content: flex-end; 
}
.search-again {
    padding: 0 1%;
    background-color: #f2f2f2;
    border-radius: 2px;
}
.instruction {
    display: none;
}
.swap-left:hover + .instruction {
    display: flex;
    padding: 5px;
    width: 90px;
    font-size: 8px;
    text-align: center;
    color: rgba(60, 60, 60, 0.66);
    background-color: rgba(247, 247, 247, 0.90);
    border: 1px solid var(--color-background);
    position: absolute;
    left: -50px;
    bottom: 15px;
}
.swap-right:hover + .instruction {
    display: flex;
    padding: 5px;
    width: 65px;
    font-size: 8px;
    text-align: center;
    color: rgba(60, 60, 60, 0.66);
    background-color: rgba(247, 247, 247, 0.90);
    border: 1px solid var(--color-background);
    position: absolute;
    right: -50px;
    bottom: 15px;
}
.expand {
    position: absolute;
    right: 2px;
    top: 2px;
}
.swap-left:hover, .swap-right:hover, .search-again:hover, .expand:hover, .download:hover {
    cursor: pointer;
}
p {
    margin: 1% 2%;
    font-size: 10px;
    font-weight: 200;
}
button {
    width: 20%;
    height: 15%;
    color: rgba(60, 60, 60, 0.66);
    background-color: #f2f2f2; 
    border-radius: 5px;
    border: none;
    cursor: pointer;
}


@media (min-width: 1024px) {
p {
    margin: 5px;
    font-size: 16px;
    font-weight: 300;
}
.parent {
    width: 70vw;
    margin: 1% 3% 15% 3%;
    padding: 2rem 4rem 2rem 4rem;
    opacity: .95;
    box-shadow: 0px 0px 1px 0px var(--color-background); 
    border-radius: 2px;
}
.header {
  margin-bottom: 2%;
}
.swap-left:hover + .instruction {
    width: 105px;
    font-size: 10px;
    left: -75px;
    bottom: 25px;
}
.swap-right:hover + .instruction {
    width: 80px;
    font-size: 10px;
    right: -50px;
    bottom: 25px;
} 
.expand:hover + .instruction {
    display: flex;
    padding: 5px;
    width: 60px;
    font-size: 10px;
    text-align: center;
    color: rgba(60, 60, 60, 0.66);
    background-color: rgba(247, 247, 247, 0.90);
    border: 1px solid var(--color-background);
    position: absolute;
    right: -50px;
    top: -25px;
} }

</style>