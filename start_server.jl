
<!-- Div para o Uploader -->
<div class="uploader-container" v-if="false">
    <q-uploader v-on:start="function(event) { handle_event(event, 'started') }" label="Upload CSV"
        hide-upload-btn="" auto-upload="" v-on:failed="function(event) { handle_event(event, 'failed') }"
        v-on:uploaded="function(event) { handle_event(event, 'uploaded') }" no-thumbnails=""
        v-on:removed="function(event) { handle_event(event, 'removed') }"
        v-on:finish="function(event) { handle_event(event, 'finished') }" :url="'/____/upload/' + channel_"
        v-on:rejected="function(event) { handle_event(event, 'rejected') }" max-files="10"
        v-on:uploading="function(event) { handle_event(event, 'uploading') }" multiple="multiple"
        :accept="upl_acceptext" v-on:added="function(event) { handle_event(event, 'added') }"
        style="max-width: 95%; width: 95%; margin: 0 auto;" :max-file-size="upl_maxsize">
    </q-uploader>
</div>

<q-table flat="" bordered="" :columns="Resulta_Table_Filter_Data.columns" v-model="Table_data"
title="Resultados" :data="Resulta_Table_Filter_Data.data" row-key="__id" v-if="visibility_result">
</q-table>

<q-table flat="" bordered="" :columns="Result_Table.columns" v-model="Table_data" title="Resultados"
:data="Result_Table.data" row-key="__id" v-if="visibility_result">
</q-table>

<q-table flat="" bordered="" :columns="Result_Table_Arvores_Duvidosas.columns" v-model="Table_data"
title="Classificação Árvores Duvidosas" :data="Result_Table_Arvores_Duvidosas.data" row-key="__id"
v-if="visibility_result">
</q-table>

<q-btn label="Return" v-if="false" v-on:click="Button_return = true" class="q-btn--night-theme">
</q-btn>
