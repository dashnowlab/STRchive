import { parse } from "yaml";
import metadataFile from "../../../CITATION.cff?raw";

/** read metadata about project itself from main citation file */
const metadata = parse(metadataFile);

export const version = metadata.version;
export const date = metadata["date-released"];
export const repo = "https://github.com/dashnowlab/STRchive";
