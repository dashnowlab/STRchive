import { useCallback, useRef, useState } from "react";

/** simple version of tanstack-query with status, error handling, de-duping */
export const useQuery = (
  /** async func that returns data */
  func,
) => {
  const [state, setState] = useState({ status: "" });

  /** keep track of latest query function run */
  const latest = useRef(Symbol());

  const query = useCallback(async () => {
    try {
      /** unique id for current query function run */
      const current = Symbol();
      latest.current = current;

      setState({ status: "loading" });

      /** run provided function to get data */
      const data = await func();

      /** if this query function run is still the latest */
      if (current === latest.current) setState({ status: "success", data });
      else console.debug("Stale query");
    } catch (error) {
      console.error(error);
      setState({ status: "error" });
    }
  }, [func]);

  return {
    /** query function */
    query,
    /** current status */
    status: state.status,
    /** current data */
    data: state.data,
  };
};
