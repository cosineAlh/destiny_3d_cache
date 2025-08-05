- Makefile: -O3 -mtune=native

- config: MonolithicStackCount: 1/2

- Bank.cpp:
  - refreshLatency
  - refreshDynamicEnergy
  - stackedDieCount
  - partitionGranularity
  - routingReadLatency
  - routingWriteLatency
  - routingResetLatency
  - routingSetLatency
  - routingRefreshLatency
  - routingReadDynamicEnergy
  - routingWriteDynamicEnergy
  - routingResetDynamicEnergy
  - routingSetDynamicEnergy
  - routingRefreshDynamicEnergy
  - routingLeakage
  - tsvArray

- BankWithHtree.cpp/BankWithoutHtree.cpp:
  - initialize:
    - params: _stackedDieCount, partitionGranularity, monolithicStackCount
  - func:
    ```cpp
    /* Calculate the physical signals that are required in routint */ 
    if (stackedDieCount > 1)
    ...

    /* TSV connections. */
    if (stackedDieCount){
        tsvArray.Initialize(tsv_type);
        tsvArray.CalculateArea();
        tsvArray.CalculateLatencyAndPower();
    }
    ```

- FunctionUnit.cpp:
  - func:
    - logical_effort

- InputParameter.cpp:
  - initialize:
    - params: min/maxStackLayer, partitionGranularity, localTsvProjection, globalTsvProjection, tsvRedundancy, monolithicStackCount

- Mat.cpp:
  - initialize:
    - params: _stackedDieCount, partitionGranularity, monolithicStackCount
  - func:
    ```cpp
    totalPredecoderOutputBits += 1 << numAddressRowPredecoderBlock1

    /* TSV connections. */
    if (stackedDieCount){
        tsvArray.Initialize(tsv_type);
        tsvArray.CalculateArea();
        tsvArray.CalculateLatencyAndPower();
    }

    refreshLatency;
    refreshDynamicEnergy;
    ```

- MemCell.cpp:
  - func:
    ```cpp
    void MemCell::ReadCellFromFile(){
      RetentionTime
      Temperature
    }

    void MemCell::ApplyPVT(){

    }
    ```

- Precharger.cpp/SenseAmp.cpp:
  - refreshLatency
  - refreshDynamicEnergy

- Result.cpp:
  - some TSV outputs

- RowDecoder.cpp:
  - func:
    ```cpp
    /* For DRAM types account for overdriven wordline. */
    readDynamicEnergy;
    ```

- SubArray.cpp:
  - initialize:
    - params: _num3DLevels
  - func:
    ```cpp
    maxWordlineCurrent;
    maxBitlineCurrent;
    ```

- Technology.cpp:
  - func:
    ```cpp
    EDRAM currents

    // Setup TSV params
    different nodes configs

    double Technology::tsv_resistance()
    double Technology::tsv_capacitance()
    double Technology::tsv_area()
    double Technology::WireTypeToTSVType()
    double Technology::SetLayerCount()
    ```

- TSV.cpp:
  - totally new file

- No Difference:
  - BasicDecoder.cpp
  - Comparator.cpp
  - formula.cpp
  - Mux.cpp
  - OutputDriver.cpp
  - Predecoder.cpp
  - Wire.cpp
