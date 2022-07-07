<script>
  import Card from './Card.vue'
  import FocusedResult from './FocusedResult.vue'
  export default {
    name: "Deck",
    data() {
      return {
        show: 0,
        fullscreen: false
      }
    },
    components: {
        Card,
        FocusedResult
    },
    props: ['count', 'dbnames', 'headers', 'results'],
    methods: {
      goback() { this.$emit("back") },
      maximize(n) {
          this.show = n;
          this.fullscreen = true;
          window.pageYOffset || document.documentElement.scrollTop || document.body.scrollTop;
      },
      minimize_screen() {
        this.fullscreen = false;
        window.pageYOffset || document.documentElement.scrollTop || document.body.scrollTop;
      }
    },
    emits: ["back"]
  }
</script>


<template>
  <div v-if="this.fullscreen == false" style="width: 100%"> 
    <h3 class="goback" @click="goback()"> ← </h3>
    <p class="instruction"> go back </p>  
  </div>

  <div v-if="this.fullscreen == false" class="board">
        
        <div class="spread">
            <Card v-for="n in count" :title="dbnames[n - 1]" :id=n v-on:fullscreen="maximize(n)"/>
        </div>
    
  </div>
  <FocusedResult v-if="this.fullscreen == true" :header="headers[this.show - 1]" :data="results[this.show - 1]"
                 v-on:back="minimize_screen()" /> 

</template>


<style scoped>


.board {
    width: 90vw;
    height: 100%;
    min-height: 90vh;
    display: flex;
    align-items: center;
    justify-content: center;
    background-color: transparent;
}
.spread {
    display: flex;
    align-items: center;
    justify-content: center;
    flex-wrap: wrap;
    height: 100%;
    width: 100%;
}
.goback {
    position: absolute; 
    background-color: transparent;
    padding: 0 1% 0.25% 1%;
    border-radius: 3px;
    left: 25px; 
    top: 10px;
    opacity: .8;
}
.goback:hover {
    cursor: pointer;
    background-color: #f2f2f2;
    opacity: 1;
}
.instruction {
    display: none;
}

@media (min-width: 1024px) {
.board {
    margin: 5vh 0;
    padding: 2% 0;
} 
.goback:hover + .instruction {
    display: flex;
    padding: 5px;
    width: 50px;
    font-size: 10px;
    text-align: center;
    color: rgba(60, 60, 60, 0.66);
    background-color: rgba(247, 247, 247, 0.90);
    box-shadow: 0px 0px 1px 0px var(--color-background); 
    position: absolute;
    left: 70px;
    top: 5px;
} }

</style>